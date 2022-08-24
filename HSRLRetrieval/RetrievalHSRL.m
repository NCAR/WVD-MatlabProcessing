% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August, 2022

function [HSRL] = RetrievalHSRL(Op,Paths,Data,Cal)
%
% Inputs: Op:    Full options structure
%         Paths: Structure containing full file paths needed for processing
%         Data:  MPD measured data, both photon counts and housekeeping
%         Cal:   Structure containing the calibration info from the
%                receiver scan and afterpulse correction
%
% Outputs: HSRL: Structure containing all calculated HSRL output
%
%% Extracting just the water vapor options for simplicity
Options = Op.HSRL;
%% Defining needed extra path information
Paths.PCASpec           = fullfile(Paths.Code,'TemperatureRetrieval','PCASpectra');
Paths.PCA.Wavelengths   = {'O2Offline'};    % Base wavelengths
Paths.PCA.Spectra       = {'RB'};   % Spectra to load
Paths.PCA.SpectraLabels = {'RayleighBr'}; % Name of spectra in code
%% Checking if temperature processing can be run
As   = {'O2Online';'O2Offline'};
Chan = {'Mol'     ;'Mol'      };
[Const,C.Mol,Data1D,S.Mol,Spectra,Surf,PossibleMol] = LoadAndPrepDataForRetrievals(As,Chan,Cal,Data,Op,Options,Paths);
Chan = {'Comb'    ;'Comb'     };
[C.Comb,~,S.Comb,PossibleComb] = IdentifyNeededInfo(Data,Cal,As,Chan,Const);
% Checking if processing is possible
if not(PossibleMol & PossibleComb)
    CWLogging('******* HSRL data not availible ******\n',Op,'Main')
    HSRL = []; return
end
clear As Chan PossibleComb PossibleMol
%% HSRL Pre-Process
% Reading Needed Data (Python HSRL and Receiver Scan)
Sp.Optics.Mol  = ReadSystemScanData(Spectra.PCA,S.Mol,Const);
Sp.Optics.Comb = ReadSystemScanData(Spectra.PCA,S.Comb,Const);
% Reorganizing the loaded structures (because they get loaded assuming only
% 2 channels are needed instead of 4)
Counts.Raw     = ReorganizeStruct(C);
Spectra.Optics = ReorganizeStruct(Sp.Optics);
clear C S Sp
% Normalizing scan magnitudes to same 
Spectra.Optics = NormalizeReceiverScan(Spectra.Optics);
%% HSRL Data Pre-Process
[Counts.Binned,BinInfo] = PreProcessLidarData(Counts.Raw,Options);
Counts.BGSub = BGSubtractLidarData(Counts.Binned,[],BinInfo,Options);
%% Determining the temperature and pressure from the weather station 
% Building an estimate of the atmosphere
HSRL.TimeStamp = Counts.BGSub.O2OfflineComb.TimeStamp;
HSRL.Range     = Counts.BGSub.O2OfflineComb.Range;
HSRL.TGuess    = (Surf.Temperature-0.0065.*Options.Range)';                % Kelvin
HSRL.PGuess    = (Surf.Pressure.*(Surf.Temperature./HSRL.TGuess').^-5.5)'; % Atmospheres
% Extracting estimates to be a self contained data structure
T = ExtractAtmoStructs(HSRL,'TGuess');
P = ExtractAtmoStructs(HSRL,'PGuess');
%% Calculating HSRL per Stillwell et al. 2020
Spectra.Rebuilt = BuildSpectra(Spectra.PCA,T,P,Data1D.Wavelength,Op);
[HSRL.Cmm,HSRL.Cmc,HSRL.Cam,HSRL.Cac] = CalculateHSRLEfficiencies(Spectra.Optics,Spectra.Rebuilt);
HSRL.Value = HSRLCalc(Counts.BGSub,HSRL);
% Calculate HSRL Data (Spuler's way)
[~,HSRL.Mask] = HSRLCalcSpuler(Counts.BGSub,HSRL.TGuess,HSRL.PGuess,Spectra,Const);
%% Blanking low altitude data
HSRL.Value(HSRL.Range<Options.BlankRange,:) = nan;
%% Calculating optical depth and backscatter coefficient
[HSRL.ABC,HSRL.OD] = BackscatterRatioToBackscatterCoefficient(HSRL.Value,HSRL.Range,HSRL.TGuess,HSRL.PGuess);
%% Smoothing
A = HSRL.Value;
A(HSRL.Mask) = nan;
HSRL.Smoothed  = WeightedSmooth(A,Options);
end

function [A] = ExtractAtmoStructs(HSRL,Type)
%
% Inputs: HSRL: HSRL data structure complete with atmospheric conditions
%         Type: Type of data structure to output (Options: TGuess, PGuess)
%
% Outputs: A:   Atmospheric data from the HSRL data structure
%
%%
A = HSRL;
A.Value = A.(Type);
A = rmfield(A,{'TGuess','PGuess'});
end

function [Out] = ReorganizeStruct(In)
%
% Input: In:   Data structure whose contents are goign to be reorganized
%              
% Output: Out: Reorganized data structure
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
for m = {'O2Offline'}
    A = S.([m{1},'Comb']).Transmission./S.([m{1},'Mol']).Transmission; 
    S.([m{1},'Mol']).Transmission = S.([m{1},'Mol']).Transmission.*(median(A));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spuler
function [BSR,Mask] = HSRLCalcSpuler(C,T,P,Sp,Const)
%
% Inputs: C:     Structure containing count arrays
%         T:     Temperature array
%         P:     Pressure array
%         Sp:    Structure containing calibration/optical Spectra
%         Const: Structure containing all universal constants
%
% Outputs: BSR:  Backscatter ratio calculated from Spuler's HSRL method
%          Mask: Data mask from low signal counts
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
% calculate the overall efficiency of the molecular channel
eta_mol     = reshape(eta_mol,1,1,length(eta_mol));  % dimension of range, time, spectral wavelength
eta_mol_all = trapz(RD.*1.2.*eta_mol,3)./trapz(RD.*1.2,3);
%% Calculate backscatter ratio
BSR = (C.O2OfflineComb.Counts./eta_comb)./(C.O2OfflineMol.Counts./eta_mol_all);
%% simple data mask to remove low count regions
threshold = 0.25;
Mask = C.O2OfflineMol.Counts< threshold;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the efficiencies from Stillwell et al. 2020 Equation 4.
function [Cmm,Cmc,Cam,Cac] = CalculateHSRLEfficiencies(Optics,Rb)
%
% Inputs: Optics: Structure containing receiver scan spectra
%         Rb:     Structure containing rebuilt Rayleigh-Brillouin spectra
%
% Outputs: Cmm:   Obs. efficicy-> molecular spectrum on molecular detector
%          Cmc:   Obs. efficicy-> molecular spectrum on combined detector
%          Cam:   Obs. efficicy-> aerosol spectrum on molecular detector
%          Cac:   Obs. efficicy-> aerosol spectrum on combined detector
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
end

% Calculating the backscatter ratio from Stillwell et al. 2020 Equation 6.
function [BSR] = HSRLCalc(C,HSRL)
%
% Inputs: C:    A strcture containing all observed photon counts
%         HSRL: A structure containing the needed channel efficienies
%
% Outputs: BSR: Backscatter ratio   [unitless]
%
%%
BSR = 1 - ((HSRL.Cmm.*C.O2OfflineComb.Counts.*C.O2OnlineMol.Counts -    ...
            HSRL.Cmc.*C.O2OfflineMol.Counts .*C.O2OnlineComb.Counts)./  ...
           (HSRL.Cam.*C.O2OfflineComb.Counts.*C.O2OnlineMol.Counts -    ...
            HSRL.Cac.*C.O2OfflineMol.Counts .*C.O2OnlineComb.Counts));
end

% This function is used to convert backscatter ratio (observed quantity) to
% an estimated backscatter ratio and an estimated aerosol optical depth
function [ABC,OD] = BackscatterRatioToBackscatterCoefficient(BSR,Range,T,P)
%
% Inputs: BSR:   Observed backscatter ratio                 [unitless]
%         Range: Range array of BSR                         [meters]
%         T:     Temperature used to calculate BSR          [Kelvin]
%         P:     Pressure used to calculate BSR             [Atmospheres]
%
% Outputs: ABC: Approximate aerosol backscatter coefficient [1/m/sr]
%          OD:  Approximate aerosol optical depth           [unitless]
%
%% Calculating molecular scattering coefficient [1/m/sr]
[BetaM,~] = RayleighBackscatterCoeff(770.1085e-9,P.*1013.25,T);
% Calculating Aerosol backscatter coefficient
ABC = BetaM.*(BSR-1);
%% Calculating optical depth
LidarRatio = 30;   % Guessing a lidar ratio (not super important)
Beta       = ABC;  % Saving a temp variable show NaNs can be ignored
Beta(isnan(Beta)) = 0; Beta(Beta<0) = 0;
OD   =  double(2.*LidarRatio.*cumsum(Beta).*(Range(2)-Range(1)));
end

% This function calculated the Rayleigh scattering efficiency as a function
% of the wavelength of probing radiation and the atmospheric parameters.
function [Beta,BetaTotal] = RayleighBackscatterCoeff (Lambda,Press,Temp)
%
% Inputs: Lambda:     Laser wavelength                   [meters]
%         Press:      Atmospheric pressure               [millibar]
%         Temp:       Atmospheric temperature            [Kelvin]
%
% Outputs: Beta:      The backscatter coefficient        [1/m/sr]
%          BetaTotal: The total scattering coefficient   [1/m]
%
%% Rayleigh Efficiency angular relationship
P = @(Theta) 0.7629.*(1+0.9324.*cosd(Theta).*cosd(Theta));
%% Calculating backscatter coeff
Beta = (2.938e-32).*(Press./Temp).*(1./(Lambda.^4.0117));   % Eq. (5.14)
%% Calculating the total scatter coefficient
BetaTotal = (Beta.*4.*pi)./P(180);                          % Eq. (5.15)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%