% Written By: Robert Stillwell
% Written For: NCAR
% 
function [Data,Retrievals,Options,Paths,RawData,RawTSData] = RunLoader(Date,System,Logging,ProcessHK,ProcessRetF,ProcessRetS)
%
% Inputs: Date:       String defining the date to run of the form YYYYMMDD
%         System:     String defining the system number to run of the form
%                     mpd_##
%         Logging:    String defining how data logging will be implimented
%                       'Full':   See everything
%                       'Skinny': See only comments from this function
%                       'None':   See nothing from this function or below
%         ProcessHK:  A boolean value: true runs housekeeping figures,
%                     false does not
%         ProcessRet: A boolean value: true runs retrievals, false does not
%                     
% Outputs: Data:      Structure containing all of the loaded and processed
%                     MPD data 
%          Options:   Structure containing user defined processing options
%          Paths:     Structure containing file path information
%          RawData:   Structure containing all of the loaded and lightly
%                     processed MPD data
%          RawTSData: Structure containing all of the loaded and lightly
%                     processed MPD time series data
%
%% Checking inputs and using default values if running as stand-alone
if nargin ~= 6
    Date        = '20210426';
    System      = 'mpd_05';
    Logging     = 'Skinny';
    ProcessHK   = false;
    ProcessRetF = false;
    ProcessRetS = false;
end
%% Adding path to utilities, defining options, and defining path info
for el = {'Definitions','HSRLRetrieval','MPDUtilities','Plotting','TemperatureRetrieval','Utilities','WVRetrieval'}
    addpath(fullfile(pwd,el{1,1}))
end; clear el
Options = DefineOptions(Date,System,Logging,ProcessHK,ProcessRetF|ProcessRetS);
[Paths,Server] = DefinePaths(Date,Options,System);
%% Reading data and pre-processing 
% Reading the calval files
CWLogging('---------Loading Cal Val files--------\n',Options,'Main')
CalInfo = ReadMPDJSonFiles(Paths.CalVal,Options.Date);
% Determining the file structure and reading the files
CWLogging('-------------Loading Data-------------\n',Options,'Main')
[RawData,AnyData] = ReadMPDData(Paths.Data,Options);
if not(AnyData)
    CWLogging('--------No data found: exiting--------\n',Options,'Main')
    Data = []; Retrievals = []; RawData = []; RawTSData = [];
    return
end
% Removing bad data
CWLogging('------------Remove bad data-----------\n',Options,'Main')
RawData = RemoveBadData(RawData);
% Force timestamps to increase monotonically
CWLogging('--Checking for monotonic time stamps--\n',Options,'Main')
RawData = CheckMonotonicTimeStamps(RawData);
% Applying wavemeter offset if possible
if isfield(CalInfo,'WMOffset') && isfield(RawData,'Laser')
    RawData.Laser = RecursiveOverwriteField(RawData.Laser,'WavelengthActual',[],-CalInfo.WMOffset);
end
% Unpacking the container/etalon/laser/MCS data to be useful
CWLogging('------------Unpack raw data-----------\n',Options,'Main')
[RawTSData, Data.Lidar] = UnpackRawData(RawData);
% Recursively looking for data breaks and marking them accordingly 
CWLogging('---------Finding data breaks----------\n',Options,'Main')
Data.TimeSeries = RecursivelyIdentifyBreaks(RawTSData,Options.BreakSize);
% Recursively pushing all 1-d data to a constant grid 
CWLogging('---------Interpolate 1d data----------\n',Options,'Main')
Data.TimeSeries = RecursivelyInterpolateStructure(Data.TimeSeries,Options.TimeGrid1d,[],Options.InterpMethod,true);
% Making sure that no time series elements are NaNs
Data.TimeSeries = RecursiveOverwriteField(Data.TimeSeries,'TimeStamp',Options.TimeGrid1d,[]);
% Loading calibration scan information 
CWLogging('----------Loading Scan file-----------\n',Options,'Main')
CalInfo.ScanData = ReadMPDCalScanFile(fieldnames(Data.Lidar.Raw),fullfile(Paths.CalFiles,CalInfo.ScanFile));
CWLogging('-------Loading Afterpulse file--------\n',Options,'Main')
CalInfo.AfterpulseData = ReadMPDAfterpulseFile(fieldnames(Data.Lidar.Raw),fullfile(Paths.CalFiles,CalInfo.AfterpulseFile));
%% Plotting field catalog infomation
if ProcessHK
    CWLogging('--------Plotting status figure--------\n',Options,'Main')
    [~,FigNum] = PlotStatusFigure(Data,RawData,Options);
    SaveUpload(FigNum,Options,Paths,Server,'Status');
    CWLogging('-----Plotting housekeeping figure-----\n',Options,'Main')
    FigNum = PlotHousekeepingFigure(Data,Options);
    SaveUpload(FigNum,Options,Paths,Server,'Housekeeping');
end
%% Process lidar data retrievals and plotting
if ProcessRetF || ProcessRetS
    % Push lidar data onto a constant grid
    CWLogging('--Push data to grid & AP Correcting---\n',Options,'Main')
    Data.Lidar.Interp = BinLidarData(Data.Lidar.Raw,CalInfo.AfterpulseData,Options.TimeGridLidar,Options.Default);
    % WV Retrieval
    if ProcessRetF
        CWLogging('--------Water Vapor Retrieval---------\n',Options,'Main')
        [Retrievals.WaterVapor] = RetrievalWV(Options,Paths,Data,CalInfo);
        FigNum = PlotWVQuicklook(Retrievals.WaterVapor,Options.Plot,Options);
        SaveUpload(FigNum,Options,Paths,Server,'Backscatter_WV')
    end
    % HSRL Retrieval
    if ProcessRetF
        CWLogging('------------HSRL Retrieval------------\n',Options,'Main')
        Retrievals.HSRL = RetrievalHSRL(Options,Paths,Data,CalInfo);
    else
        Retrievals.HSRL = [];
    end
    % Temperature Retrieval
    if ProcessRetS
        CWLogging('-----Running Temperature Retrieval----\n',Options,'Main')
        [Retrievals.Temperature,Retrievals.Python] = RetrievalTemperature(Options,Paths,Data,CalInfo,Retrievals);
        if not(isempty(Retrievals.Temperature))
            FigNum = PlotRetrievals(Retrievals,Retrievals.Python,Options,Data.TimeSeries.WeatherStation);
            SaveFigure(FigNum,Options,Paths,'Retrievals')
        end
    end
else
    Retrievals = [];
end
%% Plotting data dumps at the end of processing
% % CWLogging('---------------------Plotting data dump---------------------\n',Options,'Main')
% % PlotTSData(RawTSData,Data.TimeSeries)
% % PlotLidarData(Data.Lidar.Interp)
% % % Formatting figures
% % CWLogging('----------------------Formatting plots---------------------\n',Options,'Main')
% % FormatFigures

%% Saving quickload information
if Options.SaveQuickLoad
    CWLogging('---------Saving quickload data--------\n',Options,'Main')
    if ~exist(Paths.Quickload, 'dir')
       mkdir(Paths.Quickload)
    end
    save(fullfile(Paths.Quickload,[lower(erase(Options.System,'_')),'.',Date,'.Matlab.mat']), ...
                                          'CalInfo','Options','Retrievals')
end
end

% Function to simplify repeated calls when making quicklooks
function SaveUpload(FN,Op,Paths,Serv,Type)
%
% Inputs: FN:    The figure number to save and close
%         Op:    A structure containing all of the user defined options
%         Paths: A structure containing all the user defined file path info
%         Serv:  A boolean indicating if code is running on the server
%         Type:  A string containing the type of figure to be saved
%
%%
FTPFigure(FN,Op,Paths,Serv,Type)
SaveFigure(FN,Op,Paths,Type)
end
