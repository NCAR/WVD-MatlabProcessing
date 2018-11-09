




function [Counts,DataProducts] = RetrievalsH2O(Altitude,Counts,DataProducts,JSondeData,Map,Options,Paths,PulseInfo,SpatialAverage,SurfaceWeather,AverageRange)

%% Loading Hitran Data
HitranData = dlmread([Paths.Colormap,'/815nm_841nm_HITRAN_2008.csv'],',',[1 1 1676 8]);

%% Blanking lowest range bins
% blank lowest gates...not needed on HSRL
blank = nan.*ones(size(Counts.Integrated{1,1}(:,1:JSondeData.BlankRange/PulseInfo.BinWidth)));
Counts.CountRate{Map.Online,1}  = single(horzcat(blank, Counts.CountRate{Map.Online,1} (:,(JSondeData.BlankRange/PulseInfo.BinWidth+1):end)));
Counts.CountRate{Map.Offline,1} = single(horzcat(blank, Counts.CountRate{Map.Offline,1} (:,(JSondeData.BlankRange/PulseInfo.BinWidth+1):end)));
clear blank

%% Spectral Line Fitting
fprintf('      H2O Retrieval: Calculating Cross Section\n')
[DataProducts.Sigma{Map.Online,1}] =  ...
    SpectralLineFittingWV(Options.flag,PulseInfo.LambdaUnique{Map.Online,1},   ...
                                       PulseInfo.LambdaNearest{Map.Online,1},  ...
                                       HitranData,                             ...
                                       Counts.CountRate{Map.Online,1},        ...
                                       Altitude.RangeOriginal,                 ...
                                       SurfaceWeather.Pressure,                ...
                                       SurfaceWeather.Temperature);
[DataProducts.Sigma{Map.Offline,1}] =  ...
    SpectralLineFittingWV(Options.flag,PulseInfo.LambdaUnique{Map.Offline,1},  ...
                                       PulseInfo.LambdaNearest{Map.Offline,1}, ...
                                       HitranData,                             ...
                                       Counts.CountRate{Map.Offline,1},        ...
                                       Altitude.RangeOriginal,                 ...
                                       SurfaceWeather.Pressure,                ...
                                       SurfaceWeather.Temperature);
                                   
%% DIAL Equation to calculate Number Density and error
fprintf('      H2O Retrieval: DIAL Retrievals\n')
[DataProducts.N,DataProducts.N_Error] =  ...
    DIALEquationNarrowlySpaced(Counts.CountRate{Map.Online,1},     ...
                               Counts.Integrated{Map.Online,1},    ...
                               Counts.CountRate{Map.Offline,1},    ...
                               Counts.Integrated{Map.Offline,1},   ...
                               Counts.Background1D{Map.Online,1},  ...
                               Counts.Background1D{Map.Offline,1}, ...
                               DataProducts.Sigma{Map.Online,1},   ...
                               DataProducts.Sigma{Map.Offline,1},  ...
                               PulseInfo.BinWidth);
                           
%% Smoothing Number Density for different range zones
fprintf('      H2O Retrieval: Water Vapor Filtering\n')
if Options.flag.gradient_filter == 1
    DataProducts.N(DataProducts.N>1E18) = nan; % use this to remove high water vapor errors
end
% Calculating spatial averaging components assuming non-constant averaging
NError = []; Navg = [];
for m=1:1:size(SpatialAverage,1)
    N_error{m,1} = real(DataProducts.N_Error./sqrt(SpatialAverage(m)));      %#ok<AGROW>
    N_avg{m,1}   = movmean(DataProducts.N,SpatialAverage(m).*2,2,'omitnan'); %#ok<AGROW>
    N_avg{m,1}(isnan(DataProducts.N)) = nan;                                 %#ok<AGROW>
    % Finding the array indices to impliment non-uniform averaging 
    StartIndex = AverageRange(m);
    if m ~= 1
        StartIndex = StartIndex + 1;
    end
    if m == size(SpatialAverage,1)
        EndIndex   = size(DataProducts.N_Error,2);
    else
        EndIndex   = AverageRange(m+1);
    end
    % Recombinging non-uniform averaged data into a single contour
    NError = [NError;N_error{m,1}(:,StartIndex:EndIndex)'];                %#ok<AGROW>
    Navg   = [Navg  ;N_avg{m,1}(:,StartIndex:EndIndex)'];                  %#ok<AGROW>
end
%smooth again at the smallest resolution to avoid boundaries
DataProducts.N_Error = movmean(NError',SpatialAverage(1).*2,2,'omitnan');
DataProducts.N_Error(isnan(NError')) = nan;
DataProducts.N_avg   = movmean(Navg',SpatialAverage(1).*2,2,'omitnan');
DataProducts.N_avg(isnan(Navg')) = nan;
clear NError Navg StartIndex EndIndex N_error N_avg m

%% Mask the Number density data based on the error, correct for range center, and add WS data at lowest gate 
DataProducts.N_Masked = DataProducts.N_avg;
DataProducts.N_Masked(DataProducts.N_avg < 0) = nan; % remove non-pysical (negative) wv regions
DataProducts.N_Masked(abs(DataProducts.N_Error./DataProducts.N_avg) > 2.00) = nan; % remove high error regions
DataProducts.N_Masked(Counts.ParsedFinalGrid{1,1}./(JSondeData.MCS.bin_duration*1e-9*JSondeData.MCS.accum) > 5E6) = nan; % remove raw counts above linear count threshold (5MC/s)
% calcuate the range lag for number density (to center in range bin)
Altitude.RangeShift  = PulseInfo.BinWidth/2; %
Altitude.RangeActual = Altitude.RangeOriginal+Altitude.RangeShift; % actual range points of data
% grid to regular (75 m) gate spacing
DataProducts = RecursivelyInterpolateStructure(DataProducts,Altitude.RangeOriginal,Altitude.RangeActual,Options.InterpMethod, Options.Extrapolation);  % grid on to standard range bins
% use the weather station to fill in the bottom gates
if Options.flag.WS == 1 
    DataProducts.N_avg(:,1)    = SurfaceWeather.NumberDensity;  % gate 1, 0 meter
    DataProducts.N_Masked(:,1) = SurfaceWeather.NumberDensity;  % gate 1, 0 meter
end

% if you want to mask the data use this
if Options.flag.mask_data == 1
  DataProducts.N_avg = DataProducts.N_Masked;
end

%% Finding optical depth
% OD is - ln(I/I.o), since offline is not the same as online it needs to
% scaled by the first few good gates -- choose 300 m to 450 m
PulseInfo.ScaleOn2Off     = nanmean(Counts.CountRate{2,1}(:,300/PulseInfo.BinWidth:450/PulseInfo.BinWidth),2)./nanmean(Counts.CountRate{1,1}(:,300/PulseInfo.BinWidth:450/PulseInfo.BinWidth),2);
DataProducts.OpticalDepth = -(log(Counts.CountRate{2,1}./bsxfun(@times, Counts.CountRate{1,1}, PulseInfo.ScaleOn2Off))); % calculate column optical depth
  
end


function [Sigma] = SpectralLineFittingWV(flag, LambdaUnique, lambda, hitran, Online_Temp_Spatial_Avg, range, Surf_P, Surf_T)
%
%
%
%
%% Filling in wavelength gaps left by wavemeter communication
% this is done by assuming the nearest good value for desired wavelength is
% the same as at the bad time...note that the wavelength needs to be
% rounded to the nearest 10th of a picometer because it can be messed up in
% temporal interpolation
lambda = round(lambda,4); 
GoodMeasurements = find((isnan(lambda) == 0)==1);
lambda = interp1(GoodMeasurements,lambda(GoodMeasurements),1:1:size(lambda,1),'nearest','extrap');

%% Defining universal constants
const.m   = 18.015E-3./6.022E23; % mass of a single water molecule
const.k_B = 1.3806488e-23;       % (J/K)
const.c   = 299792458;           % (m/s) (exact)
%% Loading Hitran data
% Finding the wavenumber bounds for the Hitran linelist
WNmin = 1/(828+4)*1e7;  % 832.0 nm converted to wavenumber
WNmax = 1/(828-4)*1e7;  % 824.0 nm converted to wavenumber
%Find lines from WNmin to WNmax to calculate voigt profile
line         = double(hitran(hitran(1:size(hitran,1),1)>WNmin & ...
                             hitran(1:size(hitran,1),1)<WNmax       , 1:size(hitran,2)));
clear hitran
% Saving the hitran linelist into an easy to read data structure
Hitran.T00     = 296;          % HITRAN reference temperature [K]
Hitran.P00     = 1;            % HITRAN reference pressure [atm]
Hitran.nu0_0   = line(:,1);    % absorption line center wavenumber from HITRAN [cm^-1]
Hitran.S0      = line(:,2);    % initial linestrength from HITRAN [cm^-1/(mol*cm^-2)]   
Hitran.gammal0 = line(:,4);    % air-broadened halfwidth at T_ref and P_ref from HITRAN [cm^-1/atm]
Hitran.gamma_s = line(:,5);    % self-broadened halfwidth at T_ref and P_ref from HITRAN [cm^-1/atm]
Hitran.E       = line(:,6);    % ground state transition energy from HITRAN [cm^-1]  
Hitran.alpha   = line(:,7);    % linewidth temperature dependence factor from HITRAN
Hitran.delta   = line(:,8);    % pressure shift from HiTRAN [cm^-1 atm^-1]

%% Calculate temperature and pressure profile
if flag.WS == 1
    T0 = nanmedian(Surf_T)+273.15;
    P0 = nanmedian(Surf_P);
else
    T0 = 273+30; % surface temperature
    P0 = 0.83;
end
% Calculating a temperature profile assuming an moist adiabatic lapse rate
T = T0-0.0065.*range;     % Units of Kelvin
% Calculating a pressure profile
P = P0.*(T0./T).^-5.5;    % units of atmospheres)

%% Code to handle multiple wavelength changes during a single day
% Pre-allocating cross section array
Sigma = nan.*zeros(size(Online_Temp_Spatial_Avg));
% Looping over all unique desired wavelengths
for l=1:length(LambdaUnique)
    % Determining the desired wavelength in wavenumbers
    Hitran.nu_on = 1/(LambdaUnique(l))*1e7;
    % Looping over all altitudes
    sigma_total = zeros(1,size(Online_Temp_Spatial_Avg,2));
    for i = 1:size(Online_Temp_Spatial_Avg,2) % calculate the absorption cross section at each range
        Hitran.nu0    = Hitran.nu0_0+Hitran.delta.*(P(i)./Hitran.P00);                             % unclear if it should be Pi/P00
        Hitran.gammal = Hitran.gammal0.*(P(i)./Hitran.P00).*((Hitran.T00./T(i)).^Hitran.alpha);    %Calculate Lorentz lineweidth at P(i) and T(i)
        % revise pressure broadened halfwidth to include correction for self broadening term
        Hitran.gammad = (Hitran.nu0).*((2.0.*const.k_B.*T(i).*log(2.0))./(const.m.*const.c^2)).^(0.5);  %Calculate HWHM Doppler linewidth at T(i)
        % Term 1 in the Voigt profile
        y = (Hitran.gammal./Hitran.gammad).*((log(2.0)).^(0.5));
        % Term 2 in the Voigt profile
        x = ((Hitran.nu_on-Hitran.nu0)./Hitran.gammad).*(log(2.0)).^(0.5);
        % Setting up Voigt convolution
        t = (-(size(line,1))/2:1:size(line,1)/2-1); %set up the integration spectral step size
        t = repmat(t',1,length(x))';
        x = repmat(x,1,length(t));
        y = repmat(y,1,length(t));
        f = (exp(-t.^2.0))./(y.^2.0+(x-t).^2.0);  % combined Voigt term 1 and 2
        % Voigt integration over all of the lines at the on and offline locations
        z = trapz(t(1,:),f,2);
        integral = z;
        % Calculate linestrength at temperature T
        S = Hitran.S0.*((Hitran.T00./T(i)).^(1.5)).*exp(1.439.*Hitran.E.*((1./Hitran.T00)-(1./T(i))));
        % Calculate the Voigt profile
        K = (y(:,1)./pi).*integral;
        % Calculate the Voigt profile absorption cross section [cm^2]
        sigmav = S.*(1./Hitran.gammad).*(((log(2.0))./pi).^(0.5)).*K; %.*far_wing_on;
        % Sum contributions from all of the surrounding lines
        sigma_total(i) = sum(sigmav);
    end
    Columns2Fill = (lambda == LambdaUnique(l));
    Sigma(Columns2Fill,:) = repmat(sigma_total,sum(Columns2Fill),1);
end
end

function [N,NError] = DIALEquationNarrowlySpaced(Online,OnlineInt,Offline,OfflineInt,OnlineBG,OfflineBG,SigmaOn,SigmaOff,BinWidth)
%
%
%
%
%
%
%%
Inside = (Online.*(circshift(Offline, [0, -1])))./((circshift(Online, [0, -1])).*Offline);
del_cross = single(1./(2.*(SigmaOn-SigmaOff).*BinWidth*100));

N           =  (del_cross.*log(Inside));
N(N == inf) = nan;

% error calculation
NError = (1/2./(SigmaOn -SigmaOff )./(BinWidth*100)...
    .*sqrt((OnlineInt+OnlineBG)./OnlineInt.^2 + (circshift(OnlineInt, [0, -1])+OnlineBG)./circshift(OnlineInt, [0, -1]).^2 + ...
           (OfflineInt+OfflineBG)./OfflineInt.^2 + (circshift(OfflineInt, [0, -1])+OfflineBG)./circshift(OfflineInt, [0, -1]).^2));
end


