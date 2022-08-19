% Written By: Robert Stillwell
% Written On: February 16, 2022
% Written For: National Center for Atmospheric Research

function [Paths,Server] = DefinePaths(Date,Op,System)
%
% Inputs: Date:   String containing desired date        (Ex: '20210426')
%         Option: Structure containing user defined options
%         System: String containing system name desired (Ex: 'mpd_05')
%
% Outputs: Paths: A structure containing all of the needed filepath
%                 information to process MPD data
%
%% Defining File Path Structure
if ismac
    DataBase = '/Volumes/MPD_Data';
    CalBase  = '/Users/stillwel/Documents/StillwellResearch/Code/Instrument/eol-lidar-calvals';
    ProcBase = '/Volumes/MPD_ProcessedData';
    Server   = false;
elseif isunix
    DataBase = '/export/fog1/rsfdata/MPD';
    CalBase  = '/export/fog1/rsfdata/MPD/calibration/eol-lidar-calvals';
    ProcBase = DataBase;
    Server   = true;
end
if Op.Temp.Bootstrap
    SubFolder = 'FullProcessingBootstrap';
else
    SubFolder = 'FullProcessing';
end
Paths.CalFiles   = fullfile(CalBase,'calfiles');
Paths.CalVal     = fullfile(CalBase,'calvals',['dial',System(end),'_calvals.json']);
Paths.Code       = pwd;
Paths.Data       = fullfile(DataBase,[System,'_data'],Date(1:4),Date);
Paths.FieldCat   = '/pub/incoming/catalog/operations';
Paths.PythonData = fullfile(ProcBase,[System,'_processed_data'],'Python',...
                                [lower(erase(System,'_')),'.',Date,'.Python.nc']);
Paths.Quickload  = fullfile(ProcBase,[System,'_processed_data'],'Quickload',SubFolder);
Paths.Quicklook  = fullfile(ProcBase,[System,'_processed_data'],'Quicklook');
end