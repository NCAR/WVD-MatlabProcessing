% Written by: Scott Spuler
% Modified by: Robert Stillwell
% Modified for: National Center For Atmospheric Research
% Modification info: Downloaded: January 17, 2017
%                    
%% Setting up runtime environment
clear all; close all; clc

%% Defining processing options
Options          = DefineOptions;
Options.System   = 'DIAL01';
Options.Location = 'RELAMPAGO';
Options.Node     = 'DIAL2';   % Keep for now because hacking jsonde files

%% Defining all file paths
DatesDesired = {'181028'};

for m=1:1:size(DatesDesired,1)
    % Parsing the dates out for processing
    Date = DatesDesired{m,1};
    % Setting up the needed filepaths
    Paths = DefinePaths(Date, Options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Processing data %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read date of file
    Paths.Date = Date;
    % Loading calibration information from JSonde file
    if strcmp(Options.Node,'DIAL1')==1
        fprintf('Current JSond info for DIAL 1 is out of date.\n')
    else
        read_dial2_calvals 
        JSondeData.MCS.accum = 14000;
        JSondeData.BlankRange = 450;
        DIALAnalysis_V01_01(JSondeData, Options, Paths)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

