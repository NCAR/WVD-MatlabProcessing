% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August, 2022
%                      Added Bootstrapping December, 2023

function [HSRL,Klett,Fernald] = RetrievalHSRL(Op,Paths,Data,Cal)
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
Paths.PCA.Spectra       = {'RB'};           % Spectra to load
Paths.PCA.SpectraLabels = {'RayleighBr'};   % Spectra name in code
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
% Reorganizing loaded structures (loaded assuming 2 channels instead of 4)
Counts.Raw     = ReorganizeStruct(C);
Spectra.Optics = ReorganizeStruct(Sp.Optics);
clear C S Sp
% Normalizing scan magnitudes to same 
Spectra.Optics = NormalizeReceiverScan(Spectra.Optics);
%% HSRL Data Pre-Process
[Counts.Binned,BinInfo] = PreProcessLidarData(Counts.Raw,Options);
%% Bootstrapping
if Options.Bootstrap
    for m=1:1:Options.BootIters
        CWLogging(['   Bootstrap iteration: ',num2str(m),'\n'],Op,'Main')
        CWLogging('        Poisson Thinning & Background Subtracting\n',Op,'Sub')
        Thinned = PoissonThinLidarData(Counts.Binned,BinInfo,Options);
        % Define guess lapse rate to define atmosphere
        GuessLapse =  -(0.0065 + rand.*(0.0098-0.0065));
        % Loop over each set of thinned data
        for n=1:1:2
            CWLogging('        HSRL pre-process\n',Op,'Sub')
            Counts.BGSub = Thinned.(['PoissThined',num2str(n)]);
            % Pre-defining data structures
            T = Counts.BGSub.O2OfflineComb; P = Counts.BGSub.O2OfflineComb;
            % Define guess atmosphere
            T.Value = Surf.Temperature'+GuessLapse.*T.Range;               % Kelvin
            P.Value = Surf.Pressure'.*((Surf.Temperature'./T.Value).^ ...
                           (Const.MolMAir.*Const.G0./Const.R./GuessLapse));% Atmospheres
            % Actually doing the nuts and bolts to retrieve backscatter ratio
            HSRLTrial{m,n} = CalculateBackscatterRatio(Counts,Data1D,Op,Options,Spectra,T,P);
        end
    end
    % Adding all bootstrap averages together
    [HSRLComb,VarComb] = CalculateVariableAndVariance(HSRLTrial,{'Value','Smoothed','Temp','Press'},{'HSRL','Var'});
    % Parsing out data for returning
    HSRL              = HSRLTrial{n,1};
    HSRL.Value        = HSRLComb.Value;
    HSRL.Smoothed     = HSRLComb.Smoothed;
    HSRL.Variance     = VarComb.Value;
    HSRL.VarianceSm   = VarComb.Smoothed;
    HSRL.MaxChange      = VarComb.ValueMaxChange;
    HSRL.MaxChangeSm    = VarComb.SmoothedMaxChange;
    %HSRL.BootStrapSteps = HSRLTrial;
    HSRL.TGuess         = HSRLComb.Temp;
    HSRL.PGuess         = HSRLComb.Press;
else
    % Background subtracting photons
    CWLogging('     Background Subtracting\n',Op,'Sub')
    Counts.BGSub = BGSubtractLidarData(Counts.Binned,[],BinInfo,Options);
    % Define guess atmosphere
    T = Counts.BGSub.O2OfflineComb; P = Counts.BGSub.O2OfflineComb;              % Pre-defining data structures
    GuessLapse = -0.008;
    T.Value = Surface.Temperature'+GuessLapse.*T.Range;                    % Kelvin
    P.Value = Surface.Pressure'.*((Surface.Temperature'./T.Value).^ ...
                           (Const.MolMAir.*Const.G0./Const.R./GuessLapse));% Atmospheres
    % Calculating HSRL per Stillwell et al. 2020
    HSRL = CalculateBackscatterRatio(Counts,Data1D,Op,Options,Spectra,T,P);
    HSRL.TGuess = T.Value;
    HSRL.PGuess = P.Value;
    HSRL.Variance = HSRL.Value .*nan;
end

%% Blanking low altitude data
HSRL.Value(HSRL.Range<Options.BlankRange,:) = nan;
%% Calculating optical depth and backscatter coefficient
[HSRL.ABC,HSRL.OD] = BackscatterRatioToBackscatterCoefficient(HSRL.Value,HSRL.Range,HSRL.TGuess,HSRL.PGuess);
%% Calculating Backscatter Ratio with Klett inversion (for comparison only)
% Applying saturation correction to raw counts observed
Raw = Counts.Binned.O2OfflineComb.Counts;
Counts.Binned.O2OfflineComb.Counts = CorrectionNonparalyzable(Raw,35e-9,Data1D.MCS.O2Offline.ProfilesPerHistogram', HSRL.Range(2)-HSRL.Range(1));
% Reapplying background subtraction
Counts.BGSub = BGSubtractLidarData(Counts.Binned,[],BinInfo,Options);
%Performing the Klett inversion
SigmaR       = KlettInversion(HSRL.Range,Counts.BGSub.O2OfflineComb.Counts,Cal,1,0.01,50,770.1e-9,P.Value,T.Value);
BetaPFrenald = FernaldInversion(HSRL.Range,Counts.BGSub.O2OfflineComb.Counts,Cal,50,770.1e-9,P.Value,T.Value);
%% Normalizing non-quantitiative retrievals
[BetaM,~] = RayleighBackscatterCoeff(770.1085e-9,HSRL.PGuess.*1013.25,HSRL.TGuess);
% Saving the Klett inversion retrievals
Klett.TimeStamp = HSRL.TimeStamp;
Klett.Range     = HSRL.Range;
Klett.BSR       = SigmaR;
[Klett.ABC,Klett.OD] = BackscatterRatioToBackscatterCoefficient(Klett.BSR,HSRL.Range,HSRL.TGuess,HSRL.PGuess);
% Saving the Fernald inversion retrievals
Fernald.TimeStamp = HSRL.TimeStamp;
Fernald.Range     = HSRL.Range;
Fernald.ABC       = BetaPFrenald;
Fernald.BSR       = (Fernald.ABC + BetaM)./BetaM;
[~,Fernald.OD] = BackscatterRatioToBackscatterCoefficient(Fernald.BSR,HSRL.Range,HSRL.TGuess,HSRL.PGuess);
%% Making HSRL mask
HSRL.Mask = zeros(size(HSRL.Value));
HSRL.Mask(HSRL.Variance./HSRL.Value > 5) = 1;
HSRL.Mask(HSRL.OD > 1)                   = 1;
HSRL.Mask = DensityFiltering(HSRL.Mask,5,0.25);
HSRL.Mask = HSRL.Mask>0.5;
end

function [HSRL] = CalculateBackscatterRatio(Counts,Data1D,Op,Options,Spectra,T,P)
%% Calculating HSRL per Stillwell et al. 2020
Spectra.Rebuilt = BuildSpectra(Spectra.PCA,T,P,Data1D.Wavelength,Op);
[HSRL.Cmm,HSRL.Cmc,HSRL.Cam,HSRL.Cac] = CalculateHSRLEfficiencies(Spectra.Optics,Spectra.Rebuilt);
% Saving output data structure
HSRL.TimeStamp = Counts.BGSub.O2OfflineComb.TimeStamp;
HSRL.Range     = Counts.BGSub.O2OfflineComb.Range;
HSRL.Value     = HSRLCalc(Counts.BGSub,HSRL);
HSRL.Smoothed  = WeightedSmooth(HSRL.Value,Options);
HSRL.Temp      = T.Value;
HSRL.Press     = P.Value;
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
    S.([m{1},'Mol']).Transmission = S.([m{1},'Mol']).Transmission.*(median(A,'omitnan'));
end
end

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
Cmm = sum(A.*Rb.O2Offline.RayleighBr,3,'omitnan'); % Molecular spectrum on molecular detector
Cmc = sum(B.*Rb.O2Offline.RayleighBr,3,'omitnan'); % Molecular spectrum on combined detector
Cam = sum(A.*Rb.O2Offline.Aerosol,3,'omitnan');    % Aerosol spectrum on molecular detector
Cac = sum(B.*Rb.O2Offline.Aerosol,3,'omitnan');    % Aerosol spectrum on combined detector
end

% Calculating the backscatter ratio from Stillwell et al. 2020 Equation 6.
function [BSR] = HSRLCalc(C,HSRL)
%
% Inputs: C:    A strcture containing all observed photon counts
%         HSRL: A structure containing the needed channel efficienies
%
% Outputs: BSR: Backscatter ratio   [unitless]
%
%% Applying Savitzky-Golay filter to the molecular counts in time and range
C.O2OfflineComb.Counts(C.O2OfflineComb.Counts<10) = nan;
C.O2OfflineMol.Counts  = smoothdata2(C.O2OfflineMol.Counts,'sgolay',{3,11});

%%
BSR = 1 - ((HSRL.Cmm.*C.O2OfflineComb.Counts - HSRL.Cmc.*C.O2OfflineMol.Counts)./  ...
           (HSRL.Cam.*C.O2OfflineComb.Counts - HSRL.Cac.*C.O2OfflineMol.Counts));
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
