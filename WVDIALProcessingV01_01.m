% Written by: Scott Spuler
% Written for: National Center For Atmospheric Research
% Modified by: Robert Stillwell
%                 
%% Setting up runtime environment
% clear; close all; clc

%% Defining processing options
Options          = DefineOptions;
Options.System   = 'DIAL05';
Options.Location = 'FL1';
Options.Node     = 'DIAL2';   % Keep for now because hacking jsonde files

%% Defining all file paths
% DatesDesired = {'181026';'181027';'181028';'181029';'181030';'181031'
%                 '181101';'181102';'181103';'181104';'181105';
%                 '181106';'181107';'181108';'181109';'181110';
%                 '181111';'181112';'181113';'181114';'181115';
%                 '181116';'181117';'181118';'181119';'181120';
%                 '181121';'181122';'181123';'181124';'181125'};
% % DatesDesired = {'181125';'181126';
% %                 '181127';'181028';'181029';'181030';'181031';
% %                 '181101';'181102';'181103';'181104';'181105';
% %                 '181106';'181107';'181108';'181109';'181110';
% %                 '181111';'181112';'181113';'181114';'181115';
% %                 '181116';'181117';'181118';'181119';'181120';
% %                 '181121';'181122';'181123';'181124';'181125';
% %                 '181126';'181127';'181128';'181129';'181130';
% %                 '181201';'181202';'181203';'181204';'181205';
% %                 '181206';'181107';'181208';'181209';'181210';
% %                 '181211';'181212';'181213';'181214';'181215';
% %                 '181216';'181217'};
DatesDesired = {'190304'};

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
        Uptime(m) = DIALAnalysis_V01_02(JSondeData, Options, Paths);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

