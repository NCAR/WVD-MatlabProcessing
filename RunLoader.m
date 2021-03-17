% Written By: Robert Stillwell
% Written For: NCAR
% 

function [Data,Options,Paths,RawData,RawTSData] = RunLoader(Date,System,Logging)
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
    Date   = '20200902';
    System = 'mpd_05';
    Logging = 'Skinny';
end
%% Adding path to recursive functional utilities
addpath(fullfile(pwd,'Utilities'))
addpath(fullfile(pwd,'Plotting'))
addpath(fullfile(pwd,'HardwareDefinitions'))
%% Defining options
%%%%%%%%%%%%%%%%%%%%%%%%%% Defining user options %%%%%%%%%%%%%%%%%%%%%%%%%%
Options.BreakSize    = 15;      % Medians allowed before marking data break
Options.Date         = Date;
Options.InterpMethod = 'linear';
Options.Logging      = Logging; % 'Full', 'Skinny', 'None'
Options.SaveFigures  = true;
Options.System       = System;
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
Paths.Data      = fullfile('/export/fog1/rsfdata/MPD',[System,'_data'],Date(1:4),Date);
Paths.Quickload = fullfile('/export/fog1/rsfdata/MPD',[System,'_processed_data'],'Quickload','V04');
Paths.Quicklook = fullfile('/export/fog1/rsfdata/MPD',[System,'_processed_data'],'Quicklook');
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
%% Plotting field catalog infomation
CWLogging('--------Plotting status figure--------\n',Options,'Main')
[~,FigNum] = PlotStatusFigure(Data,RawData,Options);
SaveFigure(FigNum,Options,Paths,'Status')
CWLogging('-----Plotting housekeeping figure-----\n',Options,'Main')
FigNum = PlotHousekeepingFigure(Data,Options);
SaveFigure(FigNum,Options,Paths,'Housekeeping')
%% Process lidar data
% Push lidar data onto a constant grid
CWLogging('-----Push lidar data to known grid----\n',Options,'Main')
Data.Lidar.Interp = BinLidarData(Data.Lidar.Raw,Options.TimeGridLidar,Options.Default);

%% Plotting lidar data


%% Plotting data dumps at the end of processing
% % CWLogging('---------------------Plotting data dump---------------------\n',Options,'Main')
% % PlotTSData(RawTSData,Data.TimeSeries)
% % PlotLidarData(Data.Lidar.Interp)
% % % Formatting figures
% % CWLogging('----------------------Formatting plots---------------------\n',Options,'Main')
% % FormatFigures

%% Saving quickload information
% % % CWLogging('---------Saving quickload data--------\n',Options,'Main')
% % % cd(Paths.Quickload)
% % % save([lower(erase(Options.System,'_')),'.',Date,'.MatlabPreload.mat'],'Data','Options','RawData','RawTSData')
% % % cd(Paths.Code)

end
