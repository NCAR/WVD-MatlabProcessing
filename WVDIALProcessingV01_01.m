% Written by: Scott Spuler
% Written for: National Center For Atmospheric Research
% Modified by: Robert Stillwell
%                 
%% Setting up runtime environment
clear; close all; clc

%% Defining processing options
Options          = DefineOptions;
Options.System   = 'DIAL05';
Options.Location = 'FL1';
Options.Node     = 'DIAL2';   % Keep for now because hacking jsonde files

%% Defining all file paths
% DatesDesired = {'181101';'181102';'181103';'181104';'181105'; 
%                 '181106';'181107';'181108'};
DatesDesired = {'181115'};

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

