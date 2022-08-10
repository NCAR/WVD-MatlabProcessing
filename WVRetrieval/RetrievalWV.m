% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created March, 2022

function [WV] = RetrievalWV(Op,Paths,Data,Cal)
%
% Inputs: Op:      Full options structure
%         Options: Temperature processing specific options
%         Paths:
%         Data:
%
% Outputs:
%
%
%% Checking if temperature processing can be run
% Pulling out and loading needed data
As   = {'WVOnline';'WVOffline'};
Chan = {'';''};
[Counts.Raw,Data1D,Scan,Possible] = IdentifyNeededInfo(Data,Cal,As,Chan);
if not(Possible)
    CWLogging('***** All WV Data not availible ******\n',Op,'Main')
    WV = []; return
end
try
   Surface             = Data.TimeSeries.WeatherStation;
   Surface.TimeStamp   = Surface.TimeStamp.*60.*60;
catch
   Surface.TimeStamp   = Op.WV.TimeStamp;
   Surface.Temperature = ones(size(Surface.TimeStamp)).*15;
   Surface.Pressure    = ones(size(Surface.TimeStamp)).*0.83;
end

%% Loading in Python data as needed
Python = LoadPythonData2(Paths.PythonData);
%% Water Vapor Pre-Process
% Extra definitions
Options                 = Op.WV;
Paths.PCASpec           = fullfile(Paths.Code,'WVRetrieval','PCASpectra');
Paths.PCA.Wavelengths   = {'WVOnline';'WVOffline'};    % Base wavelengths
Paths.PCA.Spectra       = {'WV'};   % Spectra to load
Paths.PCA.SpectraLabels = {'Absorption'}; % Name of spectra in code
% Loading data needed for processing
Const       = DefineConstants;
Spectra.PCA = ReadPCASpectra(Paths,Data1D.Wavelength,Op);
% Bin lidar data to desired analysis resolution and background substracting
[Counts.Binned,BinInfo] = PreProcessLidarData(Counts.Raw,Options);
Counts.BGSub = BGSubtractLidarData(Counts.Binned,[],BinInfo,Options);
% Downsample and interpolate ancillary data to known MPD grid
Data1D  = RecursivelyInterpolate1DStructure(Data1D,Options.TimeStamp,'linear');
Surface = RecursivelyInterpolate1DStructure(Surface,Options.TimeStamp,'linear');
Surface.Temperature = Surface.Temperature + Const.C2K;
Surface.Pressure    = Surface.Pressure./Const.MBar2Atm;

%% Calculating relative backscatter
% Scale term is (range Resolution x Shots Binned x Offline Duty Cycle)
SwitchRate = Data1D.MCS.WVOffline.ProfilesPerHistogram./(Data1D.MCS.WVOnline.ProfilesPerHistogram+Data1D.MCS.WVOffline.ProfilesPerHistogram);
Scale = Options.BinRange.*Options.BinTime.*Data1D.MCS.WVOffline.ProfilesPerHistogram.* (1 - SwitchRate);
% Calculating overlap scalign factor
Overlap = interp1(Cal.Overlap.Range,Cal.Overlap.Value,Counts.Binned.WVOffline.Range,'linear','extrap');
% Calculating relative backscatter
BS = Counts.Binned.WVOffline.Counts-Counts.BGSub.WVOffline.Background;
Rb = BS.*(Counts.Binned.WVOffline.Range.^2)./Scale'./Overlap;

%% Building an estimate of the atmosphere
T = Counts.BGSub.WVOffline; P = Counts.BGSub.WVOffline;
T.Value = Surface.Temperature'-0.0065.*T.Range;                      % Kelvin
P.Value = Surface.Pressure'.*(Surface.Temperature'./T.Value).^-5.5;  % Atmospheres

%% Rebuilding needed wv spectra
Spectra.Rebuilt = BuildSpectra(Spectra.PCA,T,P,Data1D.Wavelength,Op);

%% Performing the DIAL equation
[N,NErr] =  DIALEquation(Counts.BGSub.WVOnline.Counts,                  ...
                         Counts.BGSub.WVOffline.Counts,                 ...
                         Counts.BGSub.WVOnline.Background,              ...
                         Counts.BGSub.WVOffline.Background,             ...
                         Spectra.Rebuilt.WVOnline.AbsorptionObserved,   ...
                         Spectra.Rebuilt.WVOffline.AbsorptionObserved,  ...
                         Options.BinRange);

%% Saving data structure and smoothing
WV.TimeStamp  = Counts.BGSub.WVOnline.TimeStamp;
WV.Range      = Counts.BGSub.WVOnline.Range;
WV.Value      = N.*Const.MWV.*1000;                   % molec/m^3 to g/m^3
WV.Variance   = (NErr.*Const.MWV.*1000).^2;           % molec/m^3 to g/m^3
WV.Smoothed   = WeightedSmooth(WV.Value,Options);
WV.VarianceSm = WeightedSmooth(WV.Variance,Options);

WV.RB.TimeStamp = Counts.BGSub.WVOnline.TimeStamp;
WV.RB.Range     = Counts.Binned.WVOffline.Range;
WV.RB.Value     = Rb;

WV.Python       = Python.Online;
%% Calculating data masks
% remove non-physical water vapor values
MaskNP   = WV.Smoothed > 30;
% Removing high error regions
MaskErr1  = WV.Variance > (20^2);
%MaskErr2 = abs(sqrt(WV.Variance)./WV.Smoothed2) > 5 & WV.Variance > 3^2;
MaskErr2 = 0.*MaskErr1;
MaskErr  = MaskErr1 | MaskErr2;
% Gradient filter
MaskGrad = GradientFilter(Rb, Data1D.MCS.WVOffline, BinInfo, Options);
% RB and are not the same size (in range) so downsize Gradient filter mask
MaskGrad = MaskGrad(1:size(MaskErr,1),1:size(MaskErr,2));
% Remove high count rate regions
MaskCntR = CntRate(Counts.BGSub.WVOffline.Counts, Data1D.MCS.WVOffline, BinInfo) > 2e6;
% Combining the masks
WV.Mask = MaskNP | MaskErr | MaskGrad | MaskCntR;
% Removing data with high amounts of bad data neighbors (speckle filtering)
Mask2   = DensityFiltering(WV.Mask,5,0.5);
WV.Mask = WV.Mask | Mask2;
%% Calculating smoothed profile for plotting
PreMask = WV.Value;
PreMask(WV.Mask == 1) = nan;
WV.Smoothed2  = SmoothOld2(PreMask,Options);
% %% Plotting
% PlotWVMask(WV, MaskNP, MaskErr1, MaskErr2, MaskGrad, MaskCntR, Mask2)
end

function [N,NErr] = DIALEquation(On,Off,OnBG,OffBG,SOn,SOff,Bin)
%
% Inputs: On:    Array of online photon counts, dimensions: (range,time)
%         Off:   Array of offline photon counts, dimensions: (range,time)
%         OnBG:  Array of online background counts, dimensions (1,time)
%         OffBG: Array of offline background counts, dimensions (1,time)
%         Son:   Array of absorption cross sections at online wavelength
%                of water vapor, dimensions (range,time) and units of [m^2]
%         SOff:  Array of absorption cross sectionsthe offline wavelength
%                of water vapor, dimensions (range,time) and units of [m^2]
%         Bin:   Binwidth of observations with units of [m]
%
% Output: N:     Array of water vapor concentrations, dimensions 
%                (range,time) and units of [molecules/m^3] 
%         NErr:  Array of water vapor concentrations error estimates,  
%                dimensions (range,time) and units of [molecules/m^3] 
%
%% Calculating water vapor concentration and removing known bad values
N = double(1./(2.*(SOn-SOff).*Bin)).*log((On.*Shift(Off))./(Shift(On).*Off));
N = real(N);
N(N == inf | Off < 0) = nan;
%% Calculating the non-background subtracted counts
OnT  = On+OnBG;
OffT = Off+OffBG;
%% Calculating the error estimate
NErr = 1/2./(SOn-SOff)./Bin.*sqrt(OnT./On.^2 + Shift(OnT)./Shift(On).^2 + ...
                                  OffT./Off.^2 + Shift(OffT)./Shift(Off).^2);
NErr = real(NErr);
end

function [Mask] = GradientFilter(Offline, Data, BinInfo, Options)
%
% Inputs:
%
% Outputs:
%
%% Calculating the scaling terms
% The gradient limit that we want is approximately 1000. But that limit is
% relative to a particular resolution and count rate. That 1000 limit is 
% for a laser rep rate of 7 kHz and for a range resolution of 250 ns. 
RepScale   = Data.ProfilesPerHistogram.*BinInfo.BinNum(contains(BinInfo.BinDir,'TimeStamp'))./Options.BinTime./7e3;
RangeScale = Data.RangeResolution.*BinInfo.BinNum(contains(BinInfo.BinDir,'Range'))./250;

%% Calculating the gradient and mask
[~,Grad] = gradient(Offline);
Mask = Grad > (RepScale.*RangeScale.*Options.GradFilt)';
end

function [Rate] = CntRate(Offline, Data, BinInfo)
%
% Inputs:
%
% Outputs:
%
%% Calculating the scaling terms
% Here scale is the number of shots times the bin duration. The number of
% shots is averaged above so needs to be multiplied by the number of bins
% it is averaged over and same with range resolution
Scale = Data.ProfilesPerHistogram.*BinInfo.BinNum(contains(BinInfo.BinDir,'TimeStamp')).* ...
         Data.RangeResolution.*BinInfo.BinNum(contains(BinInfo.BinDir,'Range')).*1e-9;
Rate = Offline./Scale';
end

function [Shifted] = Shift(Data)
%
% Inputs: Data: An array of data to be shifted in range with dimensions 
%               (range,time)
%
% Outputs: Shifted: An array of data that has been circularly shifted in
%                   range
%
%%
Shifted = circshift(Data,[-1,0]);
end

function [Return] = DensityFiltering(Array,Points,Allowable)
%
% Inputs: Array:     An array of data to be filtered by nan density. This
%                    array will be filtered assuming each measurement is a 
%                    column.
%         Points:    The number of points used to determine what the nan
%                    density about a particular point is.
%         Allowable: The maximum allowable density of nans about the data
%                    point of interest
%
% Outputs: Mask:     An array of the size of the input Array with a nan
%                    where density of input nans is above the allowable
%                    level
%
%% Finding the starting point
ToSkip = Points./2;
Start  = floor(ToSkip);
End    = floor(ToSkip);
%% Arrays
Mask = zeros(size(Array,1)-Points+1,size(Array,2)-Points+1); % Pre-allocating data
for n=1:1:Points
    for m=1:1:Points
        Mask = Mask + Array(m:end-Points+m,n:end-Points+n);
    end
end
Mask = Mask./Points./Points;
%% Removing points with an unallowably high density of nans
Return = Mask > Allowable;
Return = [false(Start,size(Mask,2)+Start+End); 
         [false(size(Mask,1),Start),Return,false(size(Mask,1),End)];
          false(End,size(Mask,2)+Start+End)];
end

function [Navg] = SmoothOld2(N,Options)

N(1:2,:) = nan;

AverageRange   = [1;round(1500/Options.BinRange);round(2500/Options.BinRange)];
SpatialAverage = [150; 300; 600];


% Calculating spatial averaging components assuming non-constant averaging
Navg = [];
for m=1:1:length(SpatialAverage)
    Options.SmoothRange = SpatialAverage(m);
    
    N_avg = WeightedSmooth(N,Options);
    
    % Finding the array indices to impliment non-uniform averaging 
    StartIndex = AverageRange(m);
    if m ~= 1
        StartIndex = StartIndex + 1;
    end
    if m == length(SpatialAverage)
        EndIndex   = size(N,1);
    else
        EndIndex   = AverageRange(m+1);
    end
    % Recombinging non-uniform averaged data into a single contour
    Navg   = [Navg  ;N_avg(StartIndex:EndIndex,:)];                  %#ok<AGROW>
end
Navg(isnan(N)) = nan; 

Navg = movmean(Navg,SpatialAverage(1)/Options.BinRange.*2,1,'omitnan');
Navg(isnan(N)) = nan;

end

function PlotWVMask(WV, MaskNP, MaskErr1, MaskErr2, MaskGrad, MaskCntR, Mask2)
%
% Inputs: WV:       Data structure containing all of the WV data to be
%                   passed back to the main program
%         MaskNP:   Mask of non-physical water vapor values
%         MaskErr1: Mask of high-variance water vapor values
%         MaskErr2: Mask of weird variance water vapor values
%         MaskGrad: Mask calculated from the gradient filter
%         MaskCntR: Mask calculated from the count rate
%         Mask2:    Mask calculated from the speckle filter
%
% Outputs: none
%
%% Plotting masks
figure(100)
subplot(6,1,1); pcolor(WV.TimeStamp./60./60,WV.Range./1e3,double(MaskNP));
shading flat; colorbar; caxis([0,1]); title('High Value')
subplot(6,1,2); pcolor(WV.TimeStamp./60./60,WV.Range./1e3,double(MaskErr1));
shading flat; colorbar; caxis([0,1]); title('High Variance')
subplot(6,1,3); pcolor(WV.TimeStamp./60./60,WV.Range./1e3,double(MaskErr2));
shading flat; colorbar; caxis([0,1]); title('Weird Variance')
subplot(6,1,4); pcolor(WV.TimeStamp./60./60,WV.Range./1e3,double(MaskGrad));
shading flat; colorbar; caxis([0,1]); title('Gradient Filter')
subplot(6,1,5); pcolor(WV.TimeStamp./60./60,WV.Range./1e3,double(MaskCntR));
shading flat; colorbar; caxis([0,1]); title('Count Rate')
subplot(6,1,6); pcolor(WV.TimeStamp./60./60,WV.Range./1e3,double(Mask2));
shading flat; colorbar; caxis([0,1]); title('Speckle')

end