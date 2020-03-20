% Written by: Robert Stillwell
% Written for: National Center For Atmospheric Research
%                 
%% Setting up runtime environment
clear; close all; clc

%% Selecting the systems to process
% Systems = {'DIAL01';'DIAL02';'DIAL03';'DIAL04';'DIAL05'};
Systems = {'DIAL04'};

%% Selecting the days to process
DatesDesired = {'191106'};
% DatesDesired = {'190701';'190702';'190703';'190704';'190705';
%                 '190706';'190707';'190708';'190709';'190710';
%                 '190711';'190712';'190713';'190714'};

%% Defining the processing optins
Options = DefineOptions;

%% Pre-allocating data
Uptime = zeros(size(DatesDesired,1),size(Systems,1));

% Looping over the systems
for a=1:1:size(Systems,1)
    % 
    Options.System   = Systems{a};
    Options.Location = 'FL1';
    % Looping over the days
    for m=1:1:size(DatesDesired,1)
        % Setting up the needed filepaths
        Paths = DefinePaths(DatesDesired{m,1}, Options);
        % Read date of file
        Paths.Date = DatesDesired{m,1};
        % Loading calibration information from JSonde file
        JSondeData  = ReadJSonFiles(DatesDesired{m,1},Options);
        % Processing DIAL data
        Uptime(m,a) = DIALAnalysis_V01_02(JSondeData, Options, Paths);
    end
end
