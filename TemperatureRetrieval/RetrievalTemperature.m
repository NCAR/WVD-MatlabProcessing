


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
%% Applying the data masks to and background subtracting MPD data
Data2D = RecursivelyApplyMask(Data2D);
%% Downsample and interpolate ancillary data to known MPD grid
Data1D = RecursivelyInterpolate1DStructure(Data1D,Options.TimeStamp,'linear');
Data2D = RecursivelyInterpolate2DStructure(Data2D,Options.TimeStamp,Options.Range,'linear');
%% Calculating needed spectra
GuessLapse = -0.0098;
ConstProfile = (Options.Range').*(GuessLapse);
% Data2D.NCIP.Temperature.Value = repmat(ConstProfile,1,size(Data2D.NCIP.Temperature.Value,2))+283.15;
% Data2D.NCIP.Temperature.Value = repmat(ConstProfile,1,size(Data2D.NCIP.Temperature.Value,2))+Data1D.Surface.Temperature.Value';
Spectra.Rebuilt = BuildSpectra(Spectra.PCA,Data2D.NCIP.Temperature,Data2D.NCIP.Pressure,Data1D.Wavelength,Op);
%% Calculating perturbative absorption terms
[Alpha,PerterbOrders] = PerturbativeRetrieval(Constants,Counts,Data2D,Options,Spectra,Op);
%% Converting apparent absorption to temperature 
Temperature = ConvertAlpha2Temperature(Alpha,Constants,Data1D,Data2D,Options,Data1D.Surface,Spectra,Op);
%% Smoothing temperature
Smoothing2 = repmat(gaussmf(linspace(-1,1,Options.SmoothRange/Options.BinRange)',0,0.5),1,Options.SmoothTime/Options.BinTime);
Smoothing2 = Smoothing2./sum(sum(Smoothing2));
Temperature.Smoothed = filter2(Smoothing2,Temperature.Value,'same');

%% Saving data
% % % cd('/Volumes/StillwellData01/DIAL/MPD/MatlabProcessed_Temperature/MPD05')
% % % FileName = ['MPD05_TempData_',Date,'.mat'];
% % % 
% % % save(FileName,'Alpha','Counts','Data1D','Data2D','Options','PerterbOrders','Temperature')
% % % cd(Paths.Code)


end