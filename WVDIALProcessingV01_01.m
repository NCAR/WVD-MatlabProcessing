% Written by: Scott Spuler
% Written for: National Center For Atmospheric Research
% Modified by: Robert Stillwell
%                 
%% Setting up runtime environment
clear; close all; clc

% Systems = {'DIAL01';'DIclose aAL02';'DIAL03';'DIAL04';'DIAL05'};
% Systems = {'DIAL01';'DIAL02';'DIAL04';'DIAL05'};
% Systems = {'DIAL03'};
Systems = {'DIAL05'};

for a=1:1:size(Systems,1)

%% Defining processing options
Options          = DefineOptions;
Options.System   = Systems{a};
Options.Location = 'FL1';
Options.Node     = 'DIAL2';   % Keep for now because hacking jsonde files

%% Defining all file paths
DatesDesired = {'190719'};
% DatesDesired = {'190701';'190702';'190703';'190704';'190705';
%                 '190706';'190707';'190708';'190709';'190710';
%                 '190711';'190712';'190713';'190714'};

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
        Uptime(m,a) = DIALAnalysis_V01_02(JSondeData, Options, Paths);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
