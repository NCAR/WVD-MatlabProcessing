% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created December, 2020


function [AvgTemp,Var,VarSmo,VarOld,Dt,MaxChange] = RetrievalTemperature(Op,Options,Paths,Data,PythonFile)
%
%
%
%
%
%% Temperature Pre-Process
% Extra definitions
Paths.PCASpec           = fullfile(Paths.Code,'TemperatureRetrieval','PCASpectra');
Paths.PCA.Wavelengths   = {'O2Online';'O2Offline'};    % Base wavelengths
Paths.PCA.Spectra       = {'20GHzPCA';'RB20GHzPCA'};   % Spectra to load
Paths.PCA.SpectraLabels = {'Absorption';'RayleighBr'}; % Name of spectra in code
% Pulling out and loading needed data
[Counts.Raw,Data1D]     = IdentifyNeededInfo(Data);
[Data1D.Surface,Data2D] = LoadPythonData(PythonFile);
% MPD = Data2D.MPD;
% Loading data needed for processing
Constants      = DefineConstants;
Spectra.PCA    = ReadPCASpectra(Paths);
% Reading Needed Data (Python HSRL and Receiver Scan)
Spectra.Optics = ReadSystemScanData(Spectra.PCA,Constants);                % Should load calibration scan data
% Bin lidar data to desired analysis resolution
[Counts.Binned,BinInfo] = PreProcessLidarData(Counts.Raw,Options);
% Downsample and interpolate ancillary data to known MPD grid
Data1D = RecursivelyInterpolate1DStructure(Data1D,Options.TimeStamp,'linear');
Data2D = RecursivelyInterpolate2DStructure(Data2D,Options.TimeStamp,Options.Range,'linear');

%% Bootstrapping
VarSum = zeros(length(Options.Range),length(Options.TimeStamp));
VarSumSmo = zeros(length(Options.Range),length(Options.TimeStamp));
Var    = zeros(length(Options.Range),length(Options.TimeStamp));
VarSmo    = zeros(length(Options.Range),length(Options.TimeStamp));
% Looping over bootstrap iterations
for m=1:1:Options.BootIters
    CWLogging(['Bootstrap iteration: ',num2str(m),'\n'],Op,'Main')
    CWLogging('     Poisson Thinning & Background Subtracting\n',Op,'Sub')
    Thinned = PoissonThinLidarData(Counts.Binned,BinInfo,Options);
    % Loop over each set of thinned data
    for n=1:1:2
        CWLogging('     Temperature pre-process\n',Op,'Sub')
        Counts.BGSub = Thinned.(['PoissThined',num2str(n)]);
        % Define guess atmosphere (lapse rate is a random starting variable)
        GuessLapse =  -1*(0.0065 + rand.*(0.0098-0.0065));
        ConstProfile = (Options.Range').*(GuessLapse);
        Data2D.NCIP.Temperature.Value = repmat(ConstProfile,1,size(Data2D.NCIP.Temperature.Value,2))+Data1D.Surface.Temperature.Value';
        % Build needed spectra
        Spectra.Rebuilt = BuildSpectra(Spectra.PCA,Data2D.NCIP.Temperature,Data2D.NCIP.Pressure,Data1D.Wavelength,Op);
        % Perterbative retrieval
        CWLogging('     Perterbative Retrieval\n',Op,'Sub')
        [Alpha,PerterbOrders] = PerturbativeRetrieval(Constants,Counts,Data2D,Options,Spectra,Op);
        % Convert absorption to temperature
        CWLogging('     Converting to temperature\n',Op,'Sub')
        [Temperature{n,1},Dt{m,n}] = ConvertAlpha2Temperature(Alpha,Constants,Data1D,Data2D,Options,Data1D.Surface,Spectra,Op,'Bootstrap');
        % Smooth temperature
        Smoothing = MakeSmoothingKernal(Options);
        Temperature{n,1}.Smoothed    = WeightedSmooth(Temperature{n,1}.Value,Smoothing);
    end
    % Calculate temperature average and variance from both poisson thinned profiles
    CWLogging('     Calculating variance\n',Op,'Sub')
    AvgTemp{m,1} = (Temperature{1,1}.Value    + Temperature{2,1}.Value)./2;
    AvgTemp{m,2} = (Temperature{1,1}.Smoothed + Temperature{2,1}.Smoothed)./2;
	VarOld{m,1}  = Var;
    VarOld{m,2}  = VarSmo;
    
    VarSum    = VarSum    + (Temperature{1,1}.Value    - Temperature{2,1}.Value).^2;
    VarSumSmo = VarSumSmo + (Temperature{1,1}.Smoothed - Temperature{2,1}.Smoothed).^2;
    if m >= 2
        Var    = (1./(2.*(m-1))).*VarSum;
        VarSmo = (1./(2.*(m-1))).*VarSumSmo;
    else
        Var    = (1./2).*VarSum;
        VarSmo = (1./2).*VarSumSmo;
    end
    % Calculate max change in variance for each iteration
    MaxChange(m,1) = max(max(abs(Var - VarOld{m,1})));
    MaxChange(m,2) = max(max(abs(VarSmo - VarOld{m,2})));
 
end

end