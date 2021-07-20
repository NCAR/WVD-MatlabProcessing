% Written By: Robert Stillwell
% Written For: NCAR
% 
function [Data,Retrievals,Opts,Paths,RawData,RawTSData] = RunLoader(Date,System,Logging)
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
%          Opts:      Structure containing user defined processing options
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
Opts.BreakSize     = 15;      % Medians allowed before marking databreak
Opts.Date          = Date;
Opts.InterpMethod  = 'linear';
Opts.Logging       = Logging; % 'Full', 'Skinny', 'None'
Opts.UploadFig     = true;
Opts.SaveFigures   = true;
Opts.SaveQuickLoad = false;  
Opts.System        = System;
% Temperature retrieval options
Opts.Temp.BackgroundInd = 50;     % How many pre-integration bins to   
                                  % use to estimate background noise
Opts.Temp.BinRange    = 2*37.5;   % Desired data range resolution          [meters]
Opts.Temp.BinTime     = 5*60;     % Desired data time resolution           [seconds]
Opts.Temp.SmoothRange = 300;      % Desired smoothing range res            [meters]
Opts.Temp.SmoothTime  = 30*60;    % Desired smoothing time res             [seconds]
Opts.Temp.MaxRange    = 6e3;      % Max range to run retrievals to         [meters]
Opts.Temp.MaxTime     = 24*60*60; % Max time to run retrievals to          [seconds]
Opts.Temp.MinRange    = 150;                     % Start of retrievals     [meters] 
Opts.Temp.MinTime     = Opts.Temp.BinTime./2; % Start of retrievals        [seconds]
Opts.Temp.Range       = Opts.Temp.MinRange:Opts.Temp.BinRange:Opts.Temp.MaxRange;
Opts.Temp.TimeStamp   = Opts.Temp.MinTime:Opts.Temp.BinTime:Opts.Temp.MaxTime;
%%%%%%%%%%%%%%%%%%%%%%%% Defining default options %%%%%%%%%%%%%%%%%%%%%%%%%
Opts.Default.RangeRes = 250;                    % Units are nanosceconds
Opts.Default.Range    = 16e3;                   % Units are kilometers
Opts.TimeGrid1d       = ((30:60:86400)./3600)'; % Data every 60 seconds
Opts.TimeGridLidar    = ((0:60:86400)./3600)';  % Data every 60 seconds
%%%%%%%%%%%%%%%%%%%%%%%%%% Defining data to read %%%%%%%%%%%%%%%%%%%%%%%%%%    
DataNames = {'QuantumComposer';'Container';'Etalon';'Thermocouple';
             'HumiditySensor';'Laser';'MCS';'Power';'UPS';'WeatherStation';
             'Current'};  
%% Defining filepaths
DataBase         = '/export/fog1/rsfdata/MPD';
Paths.Code       = pwd;
Paths.Data       = fullfile(DataBase,[System,'_data'],Date(1:4),Date);
Paths.PythonData = fullfile(DataBase,[System,'_processed_data'],...
                                [lower(erase(System,'_')),'.',Date(3:end),'.Python.nc']);
Paths.Quickload  = fullfile(DataBase,[System,'_processed_data'],'Quickload','TempData');
Paths.Quicklook  = fullfile(DataBase,[System,'_processed_data'],'Quicklook');
clear DataBase
%% Reading data and pre-processing 
% Determining the file structure and reading the files
CWLogging('-------------Loading Data-------------\n',Opts,'Main')
RawData = ReadMPDData(DataNames,Paths.Code,Paths.Data,Opts);
clear DataNames
% Removing bad data
CWLogging('------------Remove bad data-----------\n',Opts,'Main')
RawData = RemoveBadData(RawData);
% Force timestamps to increase monotonically
CWLogging('--Checking for monotonic time stamps--\n',Opts,'Main')
RawData = CheckMonotonicTimeStamps(RawData);
% Removing specific bad data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[IsField,~] = RecursivelyCheckIsField(RawData, {'Laser','Current'});
if IsField
    RawData.Laser.Current(RawData.Laser.Current==0) = nan;
end
clear IsField
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpacking the container/etalon/laser/MCS data to be useful
CWLogging('------------Unpack raw data-----------\n',Opts,'Main')
[RawTSData, Data.Lidar] = UnpackRawData(RawData);
% Recursively looking for data breaks and marking them accordingly 
CWLogging('---------Finding data breaks----------\n',Opts,'Main')
Data.TimeSeries = RecursivelyIdentifyBreaks(RawTSData,Opts.BreakSize);
% Recursively pushing all 1-d data to a constant grid 
CWLogging('---------Interpolate 1d data----------\n',Opts,'Main')
Data.TimeSeries = RecursivelyInterpolateStructure(Data.TimeSeries,Opts.TimeGrid1d,Opts.InterpMethod); 
% Making sure that no time series elements are NaNs
Data.TimeSeries = RecursiveOverwriteField(Data.TimeSeries,'TimeStamp',Opts.TimeGrid1d);
%% Plotting field catalog infomation
CWLogging('--------Plotting status figure--------\n',Opts,'Main')
[~,FigNum] = PlotStatusFigure(Data,RawData,Opts);
FTPFigure(FigNum,Opts,Paths,'Status')
% SaveFigure(FigNum,Options,Paths,'Status')
CWLogging('-----Plotting housekeeping figure-----\n',Opts,'Main')
FigNum = PlotHousekeepingFigure(Data,Opts);
FTPFigure(FigNum,Opts,Paths,'Housekeeping')
% SaveFigure(FigNum,Options,Paths,'Housekeeping')
%% Process lidar data
% Push lidar data onto a constant grid
CWLogging('-----Push lidar data to known grid----\n',Opts,'Main')
Data.Lidar.Interp = BinLidarData(Data.Lidar.Raw,Opts.TimeGridLidar,Opts.Default);

%% WV Retrieval 
CWLogging('--------Water Vapor Retrieval---------\n',Opts,'Main')
%% HSRL Retrieval
CWLogging('------------HSRL Retrieval------------\n',Opts,'Main')
%% Temperature Retrieval 
CWLogging('-----Running Temperature Retrieval----\n',Opts,'Main')
% % % [Retrievals.Temperature,~,Retrievals.Python] = RetrievalTemperature(Opts,Opts.Temp,Paths,Data,Paths.PythonData);
Retrievals = [];
%% Plotting lidar data
% % % FigNum = PlotRetrievals(Retrievals,Retrievals.Python,Opts,Data.TimeSeries.WeatherStation);
% % % SaveFigure(FigNum,Opts,Paths,'Retrievals')

%% Plotting data dumps at the end of processing
% % CWLogging('---------------------Plotting data dump---------------------\n',Opts,'Main')
% % PlotTSData(RawTSData,Data.TimeSeries)
% % PlotLidarData(Data.Lidar.Interp)
% % % Formatting figures
% % CWLogging('----------------------Formatting plots---------------------\n',Opts,'Main')
% % FormatFigures

%% Saving quickload information
if Opts.SaveQuickLoad
    CWLogging('---------Saving quickload data--------\n',Opts,'Main')
    cd(Paths.Quickload)
    Options = Opts;
    save([lower(erase(Opts.System,'_')),'.',Date,'.Matlab.mat'], ...
               'Data','Options','Paths','RawData','RawTSData','Retrievals')
    cd(Paths.Code)
end
end
