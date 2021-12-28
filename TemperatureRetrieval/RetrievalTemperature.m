% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created December, 2020


function [Temperature,PerterbOrders,MPD] = RetrievalTemperature(Op,Options,Paths,Data,PythonFile)
%
%
%
%
%
%
%% Extra definitions
Paths.PCASpec           = fullfile(Paths.Code,'TemperatureRetrieval','PCASpectra');
Paths.PCA.Wavelengths   = {'O2Online';'O2Offline'};    % Base wavelengths
Paths.PCA.Spectra       = {'20GHzPCA';'RB20GHzPCA'};   % Spectra to load
Paths.PCA.SpectraLabels = {'Absorption';'RayleighBr'}; % Name of spectra in code

%% Pulling out and loading needed data
[Counts.Raw,Data1D]     = IdentifyNeededInfo(Data);
[Data1D.Surface,Data2D] = LoadPythonData(PythonFile);
MPD = Data2D.MPD;
%% Loading data needed for processing
Constants      = DefineConstants;
Spectra.PCA    = ReadPCASpectra(Paths);
%% Reading Needed Data (Python HSRL and Receiver Scan)
Spectra.Optics = ReadSystemScanData(Spectra.PCA,Constants);                % Should load calibration scan data
%% Bin lidar data to desired analysis resolution and background subtracting
[Counts.Binned,Counts.BGSub] = PreProcessLidarData(Counts.Raw,Options);
% %% Applying the data masks to and background subtracting MPD data
% Data2D = RecursivelyApplyMask(Data2D);
%% Downsample and interpolate ancillary data to known MPD grid
Data1D = RecursivelyInterpolate1DStructure(Data1D,Options.TimeStamp,'linear');
Data2D = RecursivelyInterpolate2DStructure(Data2D,Options.TimeStamp,Options.Range,'linear');
%% Calculating needed spectra
GuessLapse = -0.0098;
ConstProfile = (Options.Range').*(GuessLapse);
% Data2D.NCIP.Temperature.Value = repmat(ConstProfile,1,size(Data2D.NCIP.Temperature.Value,2))+Data1D.Surface.Temperature.Value';
Spectra.Rebuilt = BuildSpectra(Spectra.PCA,Data2D.NCIP.Temperature,Data2D.NCIP.Pressure,Data1D.Wavelength,Op);
%% Calculating perturbative absorption terms
[Alpha,PerterbOrders] = PerturbativeRetrieval(Constants,Counts,Data2D,Options,Spectra,Op);
%% Converting apparent absorption to temperature 
[Temperature,Dt] = ConvertAlpha2Temperature(Alpha,Constants,Data1D,Data2D,Options,Data1D.Surface,Spectra,Op,'Ops');
Temperature.TempChange = Dt;
%% Doing more retrievals with other starting conditions for error assessment/masking
[T1,Dt1] = ConvertAlpha2Temperature(Alpha,Constants,Data1D,Data2D,Options,Data1D.Surface,Spectra,Op,'Cold');
[T2,Dt2] = ConvertAlpha2Temperature(Alpha,Constants,Data1D,Data2D,Options,Data1D.Surface,Spectra,Op,'Warm');
Temperature.IC1 = T1.Value;
Temperature.TempChangeIC1 = Dt1;
Temperature.IC2 = T2.Value;
Temperature.TempChangeIC2 = Dt2;
%% Adding extra call to do PCA temperature conversion (not yet working)
% Temperature = ConvertAlpha2TemperaturePCA(Alpha,Constants,Data1D,Data2D,Options,Data1D.Surface,Spectra,Op);
%% Masking out the lowest data points
Temperature.Absorption = Alpha;
Temperature.Value(Temperature.Range<450,:) = nan;
%% Smoothing temperature
Smoothing = MakeSmoothingKernal(Options);
Temperature.Smoothed    = WeightedSmooth(Temperature.Value,Smoothing);
Temperature.SmoothedIC1 = WeightedSmooth(Temperature.IC1,Smoothing);
Temperature.SmoothedIC2 = WeightedSmooth(Temperature.IC2,Smoothing);
end