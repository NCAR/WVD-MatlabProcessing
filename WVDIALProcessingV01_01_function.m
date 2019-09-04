

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
Options.Location = 'FL1';

%% Overwriting default options with data inputs
Options.flag.save_quicklook = save_quicklook;
Options.flag.save_data      = save_data;
Options.flag.save_netCDF    = save_netCDF;
Options.flag.save_catalog   = save_catalog;
Options.System              = node;
DatesDesired                = num2str(files);

%% Defining all file paths
Paths = DefinePaths(DatesDesired,Options);

%% Processing data 
% Read date of file
Paths.Date = DatesDesired;
% Loading calibration information from JSonde file
JSondeData = ReadJSonFiles(DatesDesired,Options);
% Processing DIAL data
DIALAnalysis_V01_02(JSondeData, Options, Paths);

end
