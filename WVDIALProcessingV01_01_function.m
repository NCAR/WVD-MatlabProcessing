

function WVDIALProcessingV01_01_function(files, save_quicklook, save_data, save_netCDF, save_catalog, node)
%
%
%
%
%
%
%% Closing all figures
close all; 

%% Defining processing options
Options          = DefineOptions;
if strcmp(node,'DIAL01')
	Options.Location = 'FL1';
else
	Options.Location = 'FL1';
end
Options.Node     = 'DIAL2';   % Keep for now because hacking jsonde files

%% Overwriting default options with data inputs
Options.flag.save_quicklook = save_quicklook;
Options.flag.save_data      = save_data;
Options.flag.save_netCDF    = save_netCDF;
Options.flag.save_catalog   = save_catalog;
Options.System              = node;
DatesDesired                = num2str(files);

%% Defining all file paths
% Parsing the dates out for processing
Date = DatesDesired;
% Setting up the needed filepaths
Paths = DefinePaths(Date,Options);
Paths.Code          = pwd; % get the current path
if strcmp(node,'DIAL01')
    Paths.Catalog       = '/pub/incoming/catalog/operations';
else
    Paths.Catalog       = '/pub/incoming/catalog/operations';
end
Paths.Figures       = ['/scr/eldora1/wvdial_',Options.System(6),'_processed_data/Quicklook'];
Paths.FigureType    = Options.System;
Paths.RawNetCDFData = ['/scr/eldora1/wvdial_',Options.System(6),'_data/20',Date(1:2),'/20',Date];

%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing data %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read date of file
Paths.Date = Date;
% Loading calibration information from JSonde file
if strcmp(Options.Node,'DIAL1')==1
    fprintf('Current JSond info for DIAL 1 is out of date.\n')
else
    read_dial2_calvals
    JSondeData.BlankRange = 450;
    DIALAnalysis_V01_02(JSondeData, Options, Paths);
end
end
