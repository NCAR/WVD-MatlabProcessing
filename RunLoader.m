% Written By: Robert Stillwell
% Written For: NCAR
% 

function [Data,Retrievals,Options,Paths,RawData,RawTSData] = RunLoader(Date,System,Logging)
%
% Inputs: Date:       String defining the date to run of the form YYYYMMDD
%         System:     String defining the system number to run of the form
%                     mpd_##
%         Logging:    String defining how data logging will be implimented
%                       'Full':   See everything
%                       'Skinny': See only comments from this function
%                       'None':   See nothing from this function or below
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
if nargin ~= 3
    Date   = '20210426';
    System = 'mpd_05';
    Logging = 'Skinny';
end
%% Adding path to recursive functional utilities
addpath(fullfile(pwd,'Utilities'))
addpath(fullfile(pwd,'Plotting'))
addpath(fullfile(pwd,'HardwareDefinitions'))
addpath(fullfile(pwd,'TemperatureRetrieval'))
%% Defining options
%%%%%%%%%%%%%%%%%%%%%%%%%% Defining user options %%%%%%%%%%%%%%%%%%%%%%%%%%
Options.BreakSize    = 15;       % Medians allowed before marking databreak
Options.Date         = Date;
Options.InterpMethod = 'linear';
Options.Logging      = Logging;  % 'Full', 'Skinny', 'None'
Options.UploadFig    = true;
Options.SaveFigures  = true;
Options.System       = System;

% Temperature retrieval options
Options.Temp.BackgroundInd = 50;     % How many pre-integration bins to   
                                     % use to estimate background noise
Options.Temp.BinRange    = 2*37.5;   % Desired data range resolution          [meters]
Options.Temp.BinTime     = 5*60;     % Desired data time resolution           [seconds]            [seconds]
Options.Temp.SmoothRange = 300;      % Desired smoothing range res            [meters]
Options.Temp.SmoothTime  = 30*60;    % Desired smoothing time res             [seconds]
Options.Temp.MaxRange    = 6e3;      % Max range to run retrievals to         [meters]
Options.Temp.MaxTime     = 24*60*60; % Max time to run retrievals to          [seconds]
Options.Temp.MinRange    = 150;                     % Start of retrievals     [meters] 
Options.Temp.MinTime     = Options.Temp.BinTime./2; % Start of retrievals     [seconds]
Options.Temp.Range       = Options.Temp.MinRange:Options.Temp.BinRange:Options.Temp.MaxRange;
Options.Temp.TimeStamp   = Options.Temp.MinTime:Options.Temp.BinTime:Options.Temp.MaxTime;

%%%%%%%%%%%%%%%%%%%%%%%% Defining default options %%%%%%%%%%%%%%%%%%%%%%%%%
Options.Default.RangeRes = 250;                    % Units are nanosceconds
Options.Default.Range    = 16e3;                   % Units are kilometers
Options.TimeGrid1d       = ((30:60:86400)./3600)'; % Data every 60 seconds
Options.TimeGridLidar    = ((0:60:86400)./3600)';  % Data every 60 seconds
%%%%%%%%%%%%%%%%%%%%%%%%%% Defining data to read %%%%%%%%%%%%%%%%%%%%%%%%%%    
DataNames = {'QuantumComposer';'Container';'Etalon';'Thermocouple';
             'HumiditySensor';'Laser';'MCS';'Power';'UPS';'WeatherStation'};  
%% Defining filepaths
Paths.Code      = pwd;
%if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu')
    Paths.Data       = fullfile('/export/fog1/rsfdata/MPD',[System,'_data'],Date(1:4),Date);
    Paths.PythonData = fullfile('/export/fog1/rsfdata/MPD',[System,'_processed_data'],...
                                [lower(erase(System,'_')),'.',Date(3:end),'.Python.nc']);
    Paths.Quickload  = fullfile('/export/fog1/rsfdata/MPD',[System,'_processed_data'],'Quickload','TempData');
    Paths.Quicklook  = fullfile('/export/fog1/rsfdata/MPD',[System,'_processed_data'],'Quicklook');
%else
%    Paths.Data       = fullfile('/Volumes/StillwellData01/DIAL/MPD/NetCDFData',[System,'_data'],Date(1:4),Date);
%    Paths.PythonData = fullfile('/Volumes/StillwellData01/DIAL/MPD/PythonProcessed',upper(erase(System,'_')),...
%                                [lower(erase(System,'_')),'.',Date(3:end),'.Python.nc']);
%    Paths.Quicklook  = fullfile('/Volumes/StillwellData01/DIAL/MPD/Quicklooks',upper(erase(System,'_')));
%    Paths.Quickload  = fullfile('/Volumes/StillwellData01/DIAL/MPD/Quickload',upper(erase(System,'_')));
%end
%% Reading data and pre-processing 
% Determining the file structure and reading the files
CWLogging('-------------Loading Data-------------\n',Options,'Main')
RawData = ReadMPDData(DataNames,Paths.Code,Paths.Data,Options);
clear DataNames
% Removing bad data
CWLogging('------------Remove bad data-----------\n',Options,'Main')
RawData = RemoveBadData(RawData);
% Force timestamps to increase monotonically
CWLogging('--Checking for monotonic time stamps--\n',Options,'Main')
RawData = CheckMonotonicTimeStamps(RawData);
% Removing specific bad data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[IsField,~] = RecursivelyCheckIsField(RawData, {'Laser','Current'});
if IsField
    RawData.Laser.Current(RawData.Laser.Current==0) = nan;
end
clear IsField
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpacking the container/etalon/laser/MCS data to be useful
CWLogging('------------Unpack raw data-----------\n',Options,'Main')
[RawTSData, Data.Lidar] = UnpackRawData(RawData);
% Recursively looking for data breaks and marking them accordingly 
CWLogging('---------Finding data breaks----------\n',Options,'Main')
Data.TimeSeries = RecursivelyIdentifyBreaks(RawTSData,Options.BreakSize);
% Recursively pushing all 1-d data to a constant grid 
CWLogging('---------Interpolate 1d data----------\n',Options,'Main')
Data.TimeSeries = RecursivelyInterpolateStructure(Data.TimeSeries,Options.TimeGrid1d,Options.InterpMethod); 
% Making sure that no time series elements are NaNs
Data.TimeSeries = RecursiveOverwriteField(Data.TimeSeries,'TimeStamp',Options.TimeGrid1d);
%% Plotting field catalog infomation
CWLogging('--------Plotting status figure--------\n',Options,'Main')
%[~,FigNum] = PlotStatusFigure(Data,RawData,Options);
%SaveFigure(FigNum,Options,Paths,'Status')
CWLogging('-----Plotting housekeeping figure-----\n',Options,'Main')
%FigNum = PlotHousekeepingFigure(Data,Options);
%SaveFigure(FigNum,Options,Paths,'Housekeeping')
%% Process lidar data
% Push lidar data onto a constant grid
CWLogging('-----Push lidar data to known grid----\n',Options,'Main')
Data.Lidar.Interp = BinLidarData(Data.Lidar.Raw,Options.TimeGridLidar,Options.Default);

%% WV Retrieval 
% CWLogging('--------Water Vapor Retrieval---------\n',Options,'Main')
%% HSRL Retrieval
% CWLogging('------------HSRL Retrieval------------\n',Options,'Main')
%% Temperature Retrieval 
CWLogging('-----Running Temperature Retrieval----\n',Options,'Main')
[Retrievals.Temperature,~,Retrievals.Python] = RetrievalTemperature(Options,Options.Temp,Paths,Data,Paths.PythonData);

%% Plotting lidar data
FigNum = PlotRetrievals(Retrievals,Retrievals.Python,Options,Data.TimeSeries.WeatherStation);
SaveFigure(FigNum,Options,Paths,'Retrievals')

%% Plotting data dumps at the end of processing
% % CWLogging('---------------------Plotting data dump---------------------\n',Options,'Main')
% % PlotTSData(RawTSData,Data.TimeSeries)
% % PlotLidarData(Data.Lidar.Interp)
% % % Formatting figures
% % CWLogging('----------------------Formatting plots---------------------\n',Options,'Main')
% % FormatFigures

%% Saving quickload information
CWLogging('---------Saving quickload data--------\n',Options,'Main')
cd(Paths.Quickload)
save([lower(erase(Options.System,'_')),'.',Date,'.Matlab.mat'],'Data','Options','Paths','RawData','RawTSData','Retrievals')
cd(Paths.Code)

end
