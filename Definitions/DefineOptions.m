% Written By: Robert Stillwell
% Written On: February 16, 2022
% Written For: National Center for Atmospheric Research

function [Op] = DefineOptions(Date,System,Logging,ProcessHK,ProcessRet)
%
% Inputs: Date:       String defining the date to run of the form YYYYMMDD
%         System:     String defining the system number to run of the form
%                     mpd_##
%         Logging:    String defining how data logging will be implimented
%                       'Full':   See everything
%                       'Skinny': See only comments from this function
%                       'None':   See nothing from this function or below
%         ProcessHK:  A boolean value: true runs housekeeping figures,
%                     false does not
%         ProcessRet: A boolean value: true runs retrievals, false does not 
%
% Outputs: Op:        A structure containing all user defined options
%                     desired for processing MPD data
%
%% Defining user options 
Op.BreakSize     = 15;          % Medians allowed before marking databreak
Op.Date          = Date;
Op.InterpMethod  = 'linear';
Op.Logging       = Logging;     % Possibilities: 'Full', 'Skinny', 'None'
Op.UploadFig     = ProcessHK;
Op.SaveFigures   = ProcessHK | ProcessRet;
Op.SaveQuickLoad = ProcessRet; 
Op.System        = System;

%% Defining default options
Op.Default.RangeRes = 250;                    % Units are nanosceconds
Op.Default.Range    = 16e3;                   % Units are kilometers
Op.TimeGrid1d       = ((30:60:86400)./3600)'; % Data every 60 seconds
Op.TimeGridLidar    = ((0:60:86400)./3600)';  % Data every 60 seconds    

%% Defining Temperature retrieval options
Op.Temp.BackgroundInd = 50;     % How many pre-integration bins to   
                                % use to estimate background noise
Op.Temp.BinRange    = 2*37.5;   % Desired data range resolution          [meters]
Op.Temp.BinTime     = 5*60;     % Desired data time resolution           [seconds]
Op.Temp.Bootstrap   = false;
Op.Temp.BootIters   = 50;       % Iterations to use when bootstraping
Op.Temp.SmoothRange = 300;      % Desired smoothing range res            [meters]
Op.Temp.SmoothTime  = 30*60;    % Desired smoothing time res             [seconds]
Op.Temp.MaxRange    = 6e3;      % Max range to run retrievals to         [meters]
Op.Temp.MaxTime     = 24*60*60; % Max time to run retrievals to          [seconds]
Op.Temp.MinRange    = 150;                     % Start of retrievals     [meters] 
Op.Temp.MinTime     = Op.Temp.BinTime./2; % Start of retrievals     [seconds]
Op.Temp.Range       = Op.Temp.MinRange:Op.Temp.BinRange:Op.Temp.MaxRange;
Op.Temp.TimeStamp   = Op.Temp.MinTime:Op.Temp.BinTime:Op.Temp.MaxTime;

%% Defining data to read 
Op.DataNames = {'Container';'Current';'Etalon';'HumiditySensor';'Laser';'MCS';
                'Power';'QuantumComposer';'Thermocouple';'UPS';'WeatherStation'}; 
end