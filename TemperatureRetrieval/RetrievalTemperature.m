% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created December, 2020

function [Temp,MPD] = RetrievalTemperature(Op,Paths,Data,Cal,Retrievals)
%
% Inputs: Op:      Full options structure
%         Options: Temperature processing specific options
%         Paths:
%         Data:
%
% Outputs:
%
%% Extracting just the temperature options for simplicity
Options = Op.Temp;
%% Defining needed extra path information
Paths.PCASpec           = fullfile(Paths.Code,'TemperatureRetrieval','PCASpectra');
Paths.HitranSpec        = fullfile(Paths.Code,'TemperatureRetrieval','HitranData','62503f11.par');
Paths.PCA.Wavelengths   = {'O2Online';'O2Offline'};    % Base wavelengths
Paths.PCA.Spectra       = {'O2';'RB'};   % Spectra to load
Paths.PCA.SpectraLabels = {'Absorption';'RayleighBr'}; % Name of spectra in code
%% Checking if temperature processing can be run
% Pulling out and loading needed data
As   = {'O2Online';'O2Offline'};
Chan = {'Comb';    'Comb'};
[Const,Counts.Raw,Data1D,Scan,Spectra,Surface,Possible] = LoadAndPrepDataForRetrievals(As,Chan,Cal,Data,Op,Options,Paths);
if not(Possible)
    CWLogging('*** Temperature data not availible ***\n',Op,'Main')
    Temp = []; MPD = []; return
end
Data1D.Surface.Temperature = BuildSimpleStruct(Surface,'Temperature');
Data1D.Surface.Pressure    = BuildSimpleStruct(Surface,'Pressure');

% Loading python data for HSRL and WV data or using onboard
if strcmp(Options.HSRLType,'Py') || strcmp(Options.HSRLType,'PyP')
    [~,Data2D,Found] = LoadPythonData(Paths.PythonData,Op);
    Data1D  = RecursivelyInterpolate1DStructure(Data1D,Options.TimeStamp,'linear');
    MPD = Data2D;
else
    Found = false;
end
% Checking if more data handling is needed
if not(Found)
    MPD = [];
    Op.Temp.HSRLType = 'On'; Options.HSRLType = 'On';
    % Mimicking NCIP information
    Data2D.Guess.Temperature   = BuildSimpleStruct(Retrievals.HSRL,'TGuess');
    Data2D.Guess.Pressure      = BuildSimpleStruct(Retrievals.HSRL,'PGuess');
end
%% Putting info in a form handy for data files and access for processing
Data2D.Onboard.HSRL = Retrievals.HSRL;
try
    Data2D.Onboard.WV   = BuildSimpleStruct(Retrievals.WaterVapor,'Smoothed2');
catch
    % If no WV processing is availible, assume WV is 0
    Data2D.Onboard.WV   = BuildSimpleStruct(Retrievals.HSRL,'Smoothed',1,0,0);
end
%% Temperature Data Pre-Process 
% Reading Needed Data (Receiver Scan)
Spectra.Optics = ReadSystemScanData(Spectra.PCA,Scan,Const);
% Bin lidar data to desired analysis resolution
[Counts.Binned,BinInfo] = PreProcessLidarData(Counts.Raw,Options);
% Downsample and interpolate ancillary data to known MPD grid
Data2D = RecursivelyInterpolate2DStructure(Data2D,Options.TimeStamp,Options.Range,'linear');

%% Bootstrapping
if Options.Bootstrap
    % Looping over bootstrap iterations
    for m=1:1:Options.BootIters
        CWLogging(['   Bootstrap iteration: ',num2str(m),'\n'],Op,'Main')
        CWLogging('        Poisson Thinning & Background Subtracting\n',Op,'Sub')
        Thinned = PoissonThinLidarData(Counts.Binned,BinInfo,Options);
        % Pre-allocating variance variables that will later be overwritten
        if m == 1
            Var = {}; VarSum = {};
        end
        % Loop over each set of thinned data
        for n=1:1:2
            CWLogging('        Temperature pre-process\n',Op,'Sub')
            Counts.BGSub = Thinned.(['PoissThined',num2str(n)]);
            % Define guess atmosphere (lapse rate is a random starting variable)
            GuessLapse =  -1*(0.0065 + rand.*(0.0098-0.0065));
            ConstProfile = (Options.Range').*(GuessLapse);
            Data2D.NCIP.Temperature.Value = ConstProfile+Data1D.Surface.Temperature.Value';
            % Actually doing the nuts and bolts to retrieve temperature
            [Alpha{m,n},~,T{m,n},Dt{m,n}] = CalculateTemperature(Const,Counts,Data1D,Data2D,Options,Spectra,Op,GuessLapse,'Bootstrap',Paths);
        end
    end
    % Adding all bootstrap averages together
    [TComb,VarComb] = CalculateTempAndVariance(T);
    % Parsing out data for returning
    Temp             = T{n,1};
    Temp.Value       = TComb.Value;
    Temp.Smoothed    = TComb.Smoothed;
    Temp.Dt          = cell2mat(reshape(Dt,size(Dt,1)*size(Dt,2),1));
    Temp.Variance    = VarComb.Value;
    Temp.VarianceSm  = VarComb.Smoothed;
    Temp.MaxChange   = VarComb.ValueMaxChange;
    Temp.MaxChangeSm = VarComb.SmoothedMaxChange;
    Temp.BootStrapSteps = T;
    Temp.Alpha       = Alpha;
else
    % Background subtracting photons
    CWLogging('     Background Subtracting\n',Op,'Sub')
    Counts.BGSub = BGSubtractLidarData(Counts.Binned,[],BinInfo,Options);
    % Actually doing the nuts and bolts to retrieve temperature
    GuessLapse = -0.008;
    [Alpha,~,Temp,Dt] = CalculateTemperature(Const,Counts,Data1D,Data2D,Options,Spectra,Op,GuessLapse,'Standard',Paths);
    % Making the output data structure
    Temp.Dt          = Dt;
    Temp.Variance    = [];
    Temp.VarianceSm  = [];
    Temp.MaxChange   = [];
    Temp.MaxChangeSm = [];
    Temp.Alpha       = Alpha;
end
end

function [Out] = BuildSimpleStruct(In,Label,TimeMult,ValMult,ValAdd)
%% Checking if the default is to preserve the structure of convert 
if nargin == 2
    TimeMult = 1; ValMult = 1; ValAdd = 0;
end
%% Creating structure
Out.TimeStamp = In.TimeStamp.*TimeMult;
if RecursivelyCheckIsField(In, {'Range'})
    Out.Range     = In.Range;
end
Out.Value     = In.(Label).*ValMult + ValAdd;
end

function [Alpha,POrders,T,Dt] = CalculateTemperature(Const,Counts,Data1D,Data2D,Options,Spectra,Op,GuessLapse,Type,Paths)
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
%% Checking which data inputs to use
if strcmp(Options.HSRLType,'Py') || strcmp(Options.HSRLType,'PyP')
    T    = Data2D.NCIP.Temperature;
    P    = Data2D.NCIP.Pressure;
    ABC  = Data2D.MPD.BSCoefficient.Value;
else
    T    = Data2D.Guess.Temperature;
    P    = Data2D.Guess.Pressure;
    ABC  = Data2D.Onboard.HSRL.ABC;
end
%% Constructing the needed spectra for processing
Spectra.Rebuilt = BuildSpectra(Spectra.PCA,T,P,Data1D.Wavelength,Op);
%% Run the Perterbative Retrieval
CWLogging('     Perterbative Retrieval\n',Op,'Sub')
[Alpha,POrders] = PerturbativeRetrieval(Const,Counts,Data2D,Options,Spectra,Op);
%% Convert absorption to temperature
CWLogging('     Converting to temperature\n',Op,'Sub')
[T,Dt] = ConvertAlpha(Alpha,Const,Data1D,Data2D,Options,Data1D.Surface,Spectra,Op,GuessLapse,Type,Paths);
%% Speckle filtering
SpeckleMask = DensityFiltering(isnan(T.Value),3,0.5);
T.Value(SpeckleMask == 1) = nan;
%% Blanking low altitude data
T.Value(T.Range<Options.BlankRange,:) = nan;
%% Removing data in excess of the allowed backscatter coefficient
% This removal is targeted at making sure cloud data is not brought into
% the temperature data. Therefore, near clouds (i.e. where smoothed data
% sees clouds) is removed.
A = WeightedSmooth(ABC > Op.Temp.BlankBSC,Options);
T.Value(A > 0.1) = nan;
%% Smooth temperature
CWLogging('     Smoothing temperature retrieval\n',Op,'Sub')
T.Smoothed = WeightedSmooth(T.Value,Options);
end
