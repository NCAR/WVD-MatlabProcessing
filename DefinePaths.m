% Written by: Robert Stillwell
% Written for: National Center For Atmospheric Research
% Modification info: Created: October 28, 2018

function [Paths] = DefinePaths(Date, Options)
%
% Inputs: 
%
% Outputs: 
%
%% Base path to change for each different platform
Paths.Base       = '/scr/eldora1';
% Paths.Base       = '/Volumes/DATAFAT32/NCAR';

%% Other paths 
Paths.Code          = '/usr/local/home/rsfdata/git/lrose-projects-eolbase/projDir/dial/MatlabV2'; % get the current path
Paths.Colormap      = [Paths.Code,'/DataFiles'];
Paths.Catalog       = '/pub/incoming/catalog/relampago';
Paths.Figures       = [Paths.Base,'/wvdial_',Options.System(6),'_processed_data/Quicklook'];
Paths.SaveData      = [Paths.Base,'/wvdial_',Options.System(6),'_processed_data/Quickload'];
Paths.FigureType    = Options.System;
Paths.RawNetCDFData = [Paths.Base,'/wvdial_',Options.System(6),'_data/20',Date(1:2),'/20',Date];

%% Adding paths to current path to run JSon file readers
addpath('/usr/local/home/rsfdata/git/lrose-projects-eolbase/projDir/dial/MatlabV2/Calibration')
addpath('/usr/local/home/rsfdata/git/lrose-projects-eolbase/projDir/dial/MatlabV2/JSon')
% addpath([pwd,'/Calibration'])
% addpath([pwd,'/JSon'])
end
