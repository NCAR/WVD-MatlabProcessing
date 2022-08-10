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
Op.UploadFig     = ProcessHK | ProcessRet;
Op.SaveFigures   = ProcessHK | ProcessRet;
Op.SaveQuickLoad = ProcessRet; 
Op.System        = System;

%% Defining default options
Op.Default.Frequency = linspace(-10,10,250);   % Units are gigahertz
Op.Default.RangeRes  = 250;                    % Units are nanosceconds
Op.Default.Range     = 16e3;                   % Units are kilometers
Op.TimeGrid1d        = ((30:60:86400)./3600)'; % Data every 60 seconds
Op.TimeGridLidar     = ((0:60:86400)./3600)';  % Data every 60 seconds

%% Defining Water Vapor retrieval options
Op.WV.BackgroundInd = 10;     % How many pre-integration bins to
                              % use to estimate background noise
Op.WV.BinRange    = 2*37.5;   % Desired data range resolution          [meters]
Op.WV.BinTime     = 4*60;     % Desired data time resolution           [seconds]
Op.WV.BlankRange  = 450;      % Altitude below which data is blanked   [meters]
Op.WV.GradFilt    = 1000;     % Relative count rate gradient to filter [ ]
Op.WV.SmoothRange = 150;      % Desired smoothing range res            [meters]
Op.WV.SmoothTime  = 4*60;     % Desired smoothing time res             [seconds]
Op.WV.MaxRange    = 6.5e3;    % Max range to run retrievals to         [meters]
Op.WV.MaxTime     = 24*60*60; % Max time to run retrievals to          [seconds]
Op.WV.MinRange    = 150;                % Start of retrievals          [meters]
Op.WV.MinTime     = Op.WV.BinTime./2;   % Start of retrievals          [seconds]
Op.WV.Range       = Op.WV.MinRange:Op.WV.BinRange:Op.WV.MaxRange;
Op.WV.TimeStamp   = Op.WV.MinTime:Op.WV.BinTime:Op.WV.MaxTime;

%% Defining HSRL retrieval options
Op.HSRL.BackgroundInd = 10;     % How many pre-integration bins to
                                % use to estimate background noise
Op.HSRL.BinRange    = 2*37.5;   % Desired data range resolution          [meters]
Op.HSRL.BinTime     = 2*60;     % Desired data time resolution           [seconds]
Op.HSRL.BlankRange  = 250;      % Altitude below which data is blanked   [meters]

Op.HSRL.SmoothRange = 2*37.5;      % Desired smoothing range res            [meters]
Op.HSRL.SmoothTime  = 2*60;     % Desired smoothing time res             [seconds]
Op.HSRL.MaxRange    = 8e3;    % Max range to run retrievals to         [meters]
Op.HSRL.MaxTime     = 24*60*60; % Max time to run retrievals to          [seconds]
Op.HSRL.MinRange    = 0;                % Start of retrievals          [meters]
Op.HSRL.MinTime     = Op.HSRL.BinTime./2;   % Start of retrievals          [seconds]
Op.HSRL.Range       = Op.HSRL.MinRange:Op.HSRL.BinRange:Op.HSRL.MaxRange;
Op.HSRL.TimeStamp   = Op.HSRL.MinTime:Op.HSRL.BinTime:Op.HSRL.MaxTime;

%% Defining Temperature retrieval options
Op.Temp.Method      = 'LBL';    % Options: 'LBL' (line-by-line) or 'PCA'
Op.Temp.BackgroundInd = 10;     % How many pre-integration bins to use to
                                % estimate background noise
Op.Temp.BinRange    = 4*37.5;   % Desired data range resolution          [meters]
Op.Temp.BinTime     = 5*60;     % Desired data time resolution           [seconds]
Op.Temp.BlankRange  = 450;      % Altitude below which data is blanked   [meters]
Op.Temp.BlankBSC    = 5e-5;     % Backscatter coefficient above which data is blanked
Op.Temp.Bootstrap   = false;
Op.Temp.BootIters   = 20;       % Iterations to use when bootstraping
Op.Temp.SmoothRange = 300;      % Desired smoothing range res            [meters]
Op.Temp.SmoothTime  = 15*60;    % Desired smoothing time res             [seconds]
Op.Temp.TempIter    = 50;       % Iterations for temperature conversion
Op.Temp.MaxRange    = 6e3;      % Max range to run retrievals to         [meters]
Op.Temp.MaxTime     = 24*60*60; % Max time to run retrievals to          [seconds]
Op.Temp.MinRange    = 150;                % Start of retrievals          [meters]
Op.Temp.MinTime     = Op.Temp.BinTime./2; % Start of retrievals          [seconds]
Op.Temp.Range       = Op.Temp.MinRange:Op.Temp.BinRange:Op.Temp.MaxRange;
Op.Temp.TimeStamp   = Op.Temp.MinTime:Op.Temp.BinTime:Op.Temp.MaxTime;

%% Defining plotting options
Op.Plot.FontSize     = 18;

Op.Plot.RB.MaxRange  = 12e3;
Op.Plot.RB.MinRange  = 0;
Op.Plot.RB.CAxis     = [1,6];

Op.Plot.WV.MaxRange  = 6e3;
Op.Plot.WV.MinRange  = 0;
Op.Plot.WV.CAxis     = [0,10];

%% Defining data to read 
Op.DataNames = {'Container';'Current';'Etalon';'HumiditySensor';'Laser';'MCS';
                'Power';'QuantumComposer';'Thermocouple';'UPS';'WeatherStation'}; 
end
