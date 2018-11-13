% Written by: Robert Stillwell
% Written for: National Center For Atmospheric Research
% This function defines all of the needed file paths as strings in a named
% structure to allow the MPD processing program to find and save data and
% processed info.
% Modification info: Created: October 28, 2018

function [Paths] = DefinePaths(Date, Options)
%
% Inputs: Date:    A string containing the desired data for processing of
%                  the form YYMMDD. The 2000 is implied. 
%         Options: Structure containing all of the user defined processing
%                   options.
%
% Outputs: Paths:  A named structure containing all of the file path
%                  information needed for finding raw data and saving
%                  processed data and figures.
%
%% Base path to change for each different platform
Paths.Base       = '/scr/eldora1';
% Paths.Base       = '/Volumes/StillwellData01/DIAL/TDIAL';
% Paths.Base       = '/Volumes/StillwellData01/DIAL/WVDIAL';

%% Other paths 
% % % Paths.Code          = '/Volumes/StillwellData01/Code/DIALAutoProcessingCode/MatlabV2';
% % % Paths.Colormap      = [Paths.Code,'/DataFiles'];
% % % % Paths.Catalog       = '/pub/incoming/catalog/operations';
% % % Paths.Figures       = [Paths.Base,'/Quicklook'];
% % % Paths.SaveData      = [Paths.Base,'/Quickload'];
% % % Paths.FigureType    = Options.System;
% % % Paths.RawNetCDFData = [Paths.Base,'/RawData/20',Date];
% % % % Paths.RawNetCDFData = [Paths.Base,'/RawData/DIAL3_NewLabview/20',Date];

Paths.Code          = '/usr/local/home/rsfdata/git/lrose-projects-eolbase/projDir/dial/Testing'; % get the current path
Paths.Colormap      = [Paths.Code,'/DataFiles'];
Paths.Catalog       = '/pub/incoming/catalog/operations';
Paths.Figures       = [Paths.Base,'/wvdial_',Options.System(6),'_processed_data/Quicklook'];
Paths.FolderType    = 'All';
Paths.SaveData      = [Paths.Base,'/wvdial_',Options.System(6),'_processed_data/Quickload'];
Paths.FigureType    = Options.System;
Paths.RawNetCDFData = [Paths.Base,'/wvdial_',Options.System(6),'_data/20',Date(1:2),'/20',Date];

%% Adding paths to current path to run JSon file readers
addpath([Paths.Code,'/Calibration'])
addpath([Paths.Code,'/JSon'])
end
