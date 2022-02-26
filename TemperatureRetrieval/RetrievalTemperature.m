% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created December, 2020


function [Temp,Variance,Dt,MaxChange,MPD] = RetrievalTemperature(Op,Paths,Data,Cal)
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
[Counts.Raw,Data1D,Scan,Possible] = IdentifyNeededInfo(Data,Cal);
if not(Possible)
    Temp = []; Variance = []; Dt = []; MaxChange = []; MPD = []; return
end
%% Temperature Pre-Process
% Extra definitions
Options                 = Op.Temp;
Paths.PCASpec           = fullfile(Paths.Code,'TemperatureRetrieval','PCASpectra');
Paths.PCA.Wavelengths   = {'O2Online';'O2Offline'};    % Base wavelengths
Paths.PCA.Spectra       = {'20GHzPCA';'RB20GHzPCA'};   % Spectra to load
Paths.PCA.SpectraLabels = {'Absorption';'RayleighBr'}; % Name of spectra in code
% Loading python data for HSRL and WV data
[Data1D.Surface,Data2D] = LoadPythonData(Paths.PythonData);
MPD = Data2D.MPD;
% Loading data needed for processing
Const       = DefineConstants;
Spectra.PCA = ReadPCASpectra(Paths);
% Reading Needed Data (Python HSRL and Receiver Scan)
Spectra.Optics = ReadSystemScanData(Spectra.PCA,Scan,Const);                % Should load calibration scan data
% Bin lidar data to desired analysis resolution
[Counts.Binned,BinInfo] = PreProcessLidarData(Counts.Raw,Options);
% Downsample and interpolate ancillary data to known MPD grid
Data1D = RecursivelyInterpolate1DStructure(Data1D,Options.TimeStamp,'linear');
Data2D = RecursivelyInterpolate2DStructure(Data2D,Options.TimeStamp,Options.Range,'linear');

%% Bootstrapping
if Options.Bootstrap
    % Looping over bootstrap iterations
    for m=1:1:Options.BootIters
        CWLogging(['Bootstrap iteration: ',num2str(m),'\n'],Op,'Main')
        CWLogging('     Poisson Thinning & Background Subtracting\n',Op,'Sub')
        Thinned = PoissonThinLidarData(Counts.Binned,BinInfo,Options);
        % Pre-allocating variance variables that will later be overwritten
        if m == 1
            Var = []; VarSum = [];
        end
        % Loop over each set of thinned data
        for n=1:1:2
            CWLogging('     Temperature pre-process\n',Op,'Sub')
            Counts.BGSub = Thinned.(['PoissThined',num2str(n)]);
            % Define guess atmosphere (lapse rate is a random starting variable)
            GuessLapse =  -1*(0.0065 + rand.*(0.0098-0.0065));
            ConstProfile = (Options.Range').*(GuessLapse);
            Data2D.NCIP.Temperature.Value = repmat(ConstProfile,1,size(Data2D.NCIP.Temperature.Value,2))+Data1D.Surface.Temperature.Value';
            % Actually doing the nuts and bolts to retrieve temperature
            [~,~,T{n,1},Dt{m,n}] = CalculateTemperature(Const,Counts,Data1D,Data2D,Options,Spectra,Op,'Bootstrap');
        end
        % Calculating variance
        [AvgTemp{m},Var,MaxChange(m,:)] = CalculateVariance(Op,m,T,Var,VarSum);
    end
    % Parsing out data for returning
    Temp     = AvgTemp(:,1);
    Variance = Var{1,1};

else
    % Background subtracting photons
    CWLogging('     Background Subtracting\n',Op,'Sub')
    Counts.BGSub = BGSubtractLidarData(Counts.Binned,[],BinInfo,Options);
    % Actually doing the nuts and bolts to retrieve temperature
    [~,~,Temp,Dt] = CalculateTemperature(Const,Counts,Data1D,Data2D,Options,Spectra,Op,'Standard');
    Variance  = [];
    MaxChange = [];
end

end


function [Alpha,POrders,T,Dt] = CalculateTemperature(Const,Counts,Data1D,Data2D,Options,Spectra,Op,Type)
%
% Inputs:
%          Type:        'Bootstrap'
%
%
% Outputs: Alpha
%          POrders     Perterbative retrieval orders
%          T           Temperature
%          Dt
%
%% Constructing the needed spectra for processing
Spectra.Rebuilt = BuildSpectra(Spectra.PCA,Data2D.NCIP.Temperature,Data2D.NCIP.Pressure,Data1D.Wavelength,Op);
%% Run the Perterbative Retrieval
CWLogging('     Perterbative Retrieval\n',Op,'Sub')
[Alpha,POrders] = PerturbativeRetrieval(Const,Counts,Data2D,Options,Spectra,Op);
%% Convert absorption to temperature
CWLogging('     Converting to temperature\n',Op,'Sub')
[T,Dt] = ConvertAlpha2Temperature(Alpha,Const,Data1D,Data2D,Options,Data1D.Surface,Spectra,Op,Type);
%% Smooth temperature
CWLogging('     Smoothing temperature retrieval\n',Op,'Sub')
T.Smoothed = WeightedSmooth(T.Value,Options);
end

function [AvgTemp,Var,MaxChange] = CalculateVariance(Op, m, Temperature, Var, VarSum)
%
%
%
%
%% Calculate temperature average and variance from poisson thinned profiles
CWLogging('     Calculating variance\n',Op,'Sub')
% Looping over the raw and smoothed data set
n = 1;
for el = {'Value','Smoothed'}
    if m == 1
        % Pre-allocating the arrays to hold current and total variance
        VarSum{n,1} = zeros(length(Options.Range),length(Options.TimeStamp));
        Var{n,1}    = zeros(length(Options.Range),length(Options.TimeStamp));
    end
    Element = el{1};
    % Calculating the average temperature from this bootstrap iteration
    AvgTemp{1,n} = (Temperature{1,1}.(Element) + Temperature{2,1}.(Element))./2;
    % Saving the old variance
    VarOld{1,n}  = Var{n,1};
    % Calculating the variance from this bootstrap iteration
    VarSum{n,1} = VarSum{n,1} + (Temperature{1,1}.(Element) - Temperature{2,1}.(Element)).^2;
    % Updating the current guess at variance
    if m >= 2
        Var{n,1} = (1./(2.*(m-1))).*VarSum{n,1};
    else
        Var{n,1} = (1./2).*VarSum{n,1};
    end
    % Calculate max change in variance for each iteration
    MaxChange(1,n) = max(max(abs(Var{n,1} - VarOld{m,n})));
    % Updating the loop counter
    n = n + 1;
end
end

