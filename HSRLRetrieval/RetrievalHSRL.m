% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August, 2022

function [HSRL] = RetrievalHSRL(Op,Paths,Data,Cal)
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
As   = {'O2Online';'O2Offline'};
Chan = {'Mol'     ;'Mol'      };
[C.Mol ,  ~   ,S.Mol ,PossibleMol ] = IdentifyNeededInfo(Data,Cal,As,Chan);
Chan = {'Comb'    ;'Comb'     };
[C.Comb,Data1D,S.Comb,PossibleComb] = IdentifyNeededInfo(Data,Cal,As,Chan);
% Loading python data for HSRL and WV data
if not(PossibleMol & PossibleComb)
    CWLogging('**** All HSRL Data not availible *****\n',Op,'Main')
    HSRL = []; return
end
clear As Chan PossibleComb PossibleMol

%% HSRL Pre-Process
% Extra definitions
Options                 = Op.HSRL;
Paths.PCASpec           = fullfile(Paths.Code,'TemperatureRetrieval','PCASpectra');
Paths.PCA.Wavelengths   = {'O2Online';'O2Offline'};    % Base wavelengths
Paths.PCA.Spectra       = {'RB'};   % Spectra to load
Paths.PCA.SpectraLabels = {'RayleighBr'}; % Name of spectra in code
% Loading data needed for processing
Const       = DefineConstants;
Spectra.PCA = ReadPCASpectra(Paths,Data1D.Wavelength,Op);
% Reading Needed Data (Python HSRL and Receiver Scan)
Sp.Optics.Mol  = ReadSystemScanData(Spectra.PCA,S.Mol,Const);              % Should load calibration scan data
Sp.Optics.Comb = ReadSystemScanData(Spectra.PCA,S.Comb,Const);
% Reorganizing the loaded structures
Counts.Raw     = ReorganizeStruct(C);
Spectra.Optics = ReorganizeStruct(Sp.Optics);
clear C S Sp
% Normalizing scan magnitudes to same 
Spectra.Optics = NormalizeReceiverScan(Spectra.Optics);
% Bin lidar data to desired analysis resolution
[Counts.Binned,BinInfo] = PreProcessLidarData(Counts.Raw,Options);
% Downsample and interpolate ancillary data to known MPD grid
Data1D = RecursivelyInterpolate1DStructure(Data1D,Options.TimeStamp,'linear');
% Background subtracting photons
CWLogging('     Background Subtracting\n',Op,'Sub')
Counts.BGSub = BGSubtractLidarData(Counts.Binned,[],BinInfo,Options);

%% Determining the temperature and pressure from the weather station 
try
   Surf             = Data.TimeSeries.WeatherStation;
   Surf.TimeStamp   = Surf.TimeStamp.*60.*60;
catch
   Surf.TimeStamp   = Op.WV.TimeStamp;
   Surf.Temperature = ones(size(Surf.TimeStamp)).*15;
   Surf.Pressure    = ones(size(Surf.TimeStamp)).*0.83;
end
Surf = RecursivelyInterpolate1DStructure(Surf,Options.TimeStamp,'linear');
Surf.Temperature = Surf.Temperature + Const.C2K;
Surf.Pressure    = Surf.Pressure./Const.MBar2Atm;

% Building an estimate of the atmosphere
HSRL.TimeStamp = Counts.BGSub.O2OfflineComb.TimeStamp;
HSRL.Range     = Counts.BGSub.O2OfflineComb.Range;
HSRL.TGuess    = (Surf.Temperature-0.0065.*Options.Range)';           % Kelvin
HSRL.PGuess    = (Surf.Pressure.*(Surf.Temperature./HSRL.TGuess').^-5.5)'; % Atmospheres

T = ExtractAtmoStructs(HSRL,'TGuess');
P = ExtractAtmoStructs(HSRL,'PGuess');


%% Calculate HSRL Data (Spuler's way)
[HSRL.Value,HSRL.Mask] = HSRLCalc(Counts.BGSub,HSRL.TGuess,HSRL.PGuess,Spectra,Const);

%% Calculating HSRL per Stillwell et al. 2020
Spectra.Rebuilt = BuildSpectra(Spectra.PCA,T,P,Data1D.Wavelength,Op);
[HSRL.Cmm,HSRL.Cmc,HSRL.Cam,HSRL.Cac] = CalculateHSRLEfficiencies(Spectra.Optics,Spectra.Rebuilt);
HSRL.ValueV2 = HSRLCalcV2(Counts.BGSub,HSRL);

%% Blanking low altitude data
HSRL.Value(HSRL.Range<Options.BlankRange,:) = nan;
HSRL.ValueV2(HSRL.Range<Options.BlankRange,:) = nan;
%% Smoothing
A = HSRL.Value;
A(HSRL.Mask) = nan;
HSRL.Smoothed  = WeightedSmooth(A,Options);

A = HSRL.ValueV2;
A(HSRL.Mask) = nan;
HSRL.SmoothedV2  = WeightedSmooth(A,Options);
end

function [A] = ExtractAtmoStructs(HSRL,Type)
%
%
%
%
%
%%
A = HSRL;
A.Value = A.(Type);
A = rmfield(A,{'TGuess','PGuess'});
end

function [Out] = ReorganizeStruct(In)
%
% Input: In:   
%              
% Output: Out: 
%
%%
for DType = fields(In)'
    DetectorType = DType{1};
    for LType = fields(In.(DetectorType))'
        LaserType = LType{1};
        Out.([LaserType,DetectorType]) = In.(DetectorType).(LaserType);
    end
end
end

function [S] = NormalizeReceiverScan(S)
%
% Inputs: S:  Structure containing receiver scan information
%             
% Outputs: S: Structure containing normalized receiver scan information
%
%% Normalizing scan parameters
for m = {'O2Online','O2Offline'}
    A = S.([m{1},'Comb']).Transmission./S.([m{1},'Mol']).Transmission; 
    S.([m{1},'Mol']).Transmission = S.([m{1},'Mol']).Transmission.*(median(A));
end
end


% Spuler
function [BSR,Mask] = HSRLCalc(C,T,P,Sp,Const)
%
% Inputs: C:     Structure containing count arrays
%         T:     Temperature array
%         P:     Pressure array
%         Sp:    Structure containing calibration/optical Spectra
%         Const: Structure containing all universal constants
%
%
%% Constants (needs to come through data but for now okay as constant)
lambda  = 770.1085*1E-9; % wavelength in nm (should really get this elsewhere)

%% Calculating the geometric overlap function
O2_geometric_correction = mean(C.O2OnlineComb.Counts./C.O2OnlineMol.Counts,2,'omitnan')';
% trying to make this an automatic calculated value, but it blows up in clouds
O2_geometric_correction(ismember(O2_geometric_correction,[0,Inf,-Inf])) = nan;

%% Calculate optical efficiencies 
try
    % the geometric correction is ignored for now
    eta_comb = median(O2_geometric_correction,'omitnan');
catch
    warning('Problem with the geometric correction, assigning a fixed value.');
    eta_comb = 0.5988;  % override, and use this value for now
end

% only a fraction of total molecular gets through the K filter
% need to know filter width, molecular scatter width RB which is range dependent
eta_comb_cal = 1./(mean(Sp.Optics.O2OfflineMol.Transmission./ ...
                        Sp.Optics.O2OfflineComb.Transmission,'omitnan'));

eta_mol = Sp.Optics.O2OfflineMol.Transmission./ ...
         (Sp.Optics.O2OfflineComb.Transmission/eta_comb_cal);

%% Calculate the Rayleigh-Doppler broadening for each height
% from Fiocco and DeWolf 1968
K0  = 2*pi/(lambda)/Const.C;
K   = K0 + K0; % backscatter is double Doppler shifted
lam = Sp.Optics.O2OfflineMol.Lambda.*1e-9;
lam = reshape(lam,1,1,length(lam));  % dimension of range, time, spectral wavelength
RD  = sqrt(Const.MAir./(2*pi*K^2*Const.Kb.*T))...
    .*exp((-1*Const.MAir./(2*K^2*Const.Kb.*T)).*(2*pi*((1./lam)-(1./lambda))).^2);
% Bosenberg 1998 suggests a simple Brillouin broadening correction: RD width x 1.2 ~ RDB width
% (valid in the lower 10km of a standard atmosphere)
RD = RD.*1.2;

%%%%%%%%%%%%%%%%%%%%% Looks like 1.2 is applied twice %%%%%%%%%%%%%%%%%%%%%

% calculate the overall efficiency of the molecular channel
eta_mol     = reshape(eta_mol,1,1,length(eta_mol));  % dimension of range, time, spectral wavelength
eta_mol_all = trapz(RD.*1.2.*eta_mol,3)./trapz(RD.*1.2,3);

%% Calculate backscatter ratio
BSR = (C.O2OfflineComb.Counts./eta_comb)./(C.O2OfflineMol.Counts./eta_mol_all);
%% simple data mask to remove low count regions
threshold = 0.25;
Mask = C.O2OfflineMol.Counts< threshold;
% BSR(C.O2OfflineMol.Counts< threshold)=nan;  % remove very low count rate molecular data

end


% Stillwell
function [Cmm,Cmc,Cam,Cac] = CalculateHSRLEfficiencies(Optics,Rb)
%
%
%
%
%% Normalizing spectra based on max value of combined receiver scan
Normalization = max(Optics.O2OfflineComb.Transmission);
Optics.O2OfflineComb.Transmission = Optics.O2OfflineComb.Transmission./Normalization;
Optics.O2OfflineMol.Transmission  = Optics.O2OfflineMol.Transmission./Normalization;
% Making transmission arrays the correct for the optics
A = permute(Optics.O2OfflineMol.Transmission,[3,1,2]);
B = permute(Optics.O2OfflineComb.Transmission,[3,1,2]);

%% Calculating the spectrally resolved efficiency of each detector
Cmm = sum(A.*Rb.O2Offline.RayleighBr,3); % Molecular spectrum on molecular detector
Cmc = sum(B.*Rb.O2Offline.RayleighBr,3); % Molecular spectrum on combined detector
Cam = sum(A.*Rb.O2Offline.Aerosol,3);    % Aerosol spectrum on molecular detector
Cac = sum(B.*Rb.O2Offline.Aerosol,3);    % Aerosol spectrum on combined detector

% % % % % Cmm = A.*Rb.O2Offline.RayleighBr; % Molecular spectrum on molecular detector
% % % % % Cmc = B.*Rb.O2Offline.RayleighBr; % Molecular spectrum on combined detector
% % % % % Cam = A.*Rb.O2Offline.Aerosol;    % Aerosol spectrum on molecular detector
% % % % % Cac = B.*Rb.O2Offline.Aerosol;    % Aerosol spectrum on combined detector
% % % % % Cmm_sum = sum(Cmm,3);
% % % % % Cmc_sum = sum(Cmc,3);
% % % % % Cam_sum = sum(Cam,3);
% % % % % Cac_sum = sum(Cac,3);
% % % % % figure(1);
% % % % % set(gcf,'Color',[1 1 1],'Position',[106,97,967,920]);
% % % % % subplot(3,1,1)
% % % % % plot(Optics.O2OfflineMol.Lambda,squeeze(Rb.O2Offline.RayleighBr(10,10,:)),'k','linewidth',2)
% % % % % xlabel('Wavelength [nm]'); title('Normalized Rayleigh-Brillouin Spectrum')
% % % % % set(gca,'xgrid','on','ygrid','on','box','on','linewidth',2,'FontSize',14)
% % % % % subplot(3,1,2); hold on;
% % % % % plot(Optics.O2OfflineMol.Lambda,squeeze(A),'k','linewidth',2); 
% % % % % plot(Optics.O2OfflineMol.Lambda,squeeze(B),'r','linewidth',2);
% % % % % xlabel('Wavelength [nm]'); title('Normalized Optical Efficiency')
% % % % % legend('Molecular Channel','Combined Channel')
% % % % % set(gca,'xgrid','on','ygrid','on','box','on','linewidth',2,'FontSize',14)
% % % % % subplot(3,1,3); hold on;
% % % % % plot(Optics.O2OfflineMol.Lambda,squeeze(Cmm(10,10,:)),'k','linewidth',2); 
% % % % % plot(Optics.O2OfflineMol.Lambda,squeeze(Cmc(10,10,:)),'r','linewidth',2);
% % % % % xlabel('Wavelength [nm]'); title('Spectrum Transmission Efficiency')
% % % % % legend('Molecular Channel','Combined Channel')
% % % % % set(gca,'xgrid','on','ygrid','on','box','on','linewidth',2,'FontSize',14)
end

function [BSR] = HSRLCalcV2(C,HSRL)
%
%
%
%
%%
BSR = 1 - ((HSRL.Cmm.*C.O2OfflineComb.Counts.*C.O2OnlineMol.Counts -    ...
            HSRL.Cmc.*C.O2OfflineMol.Counts .*C.O2OnlineComb.Counts)./ ...
           (HSRL.Cam.*C.O2OfflineComb.Counts.*C.O2OnlineMol.Counts -    ...
            HSRL.Cac.*C.O2OfflineMol.Counts .*C.O2OnlineComb.Counts));
end