% Written By: Robert Stillwell
% Written On: February 16, 2022
% Written For: National Center for Atmospheric Research

function [Op] = DefineOptions(Date,System,Logging,ProcessHK,ProcessRet,BSOverwrite)
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
%         BSOverwrite:A boolean value: true forces bootstrap processing,
%                     false allows default settings
%
% Outputs: Op:        A structure containing all user defined options
%                     desired for processing MPD data
%
%% Defining user options 
Op.Afterpulse    = false;
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
% Default settings (used to simplify definitions with common settings)
Def.BackgroundLow = 16e3;     % Altitude bottom to calcluate background [m]
Def.BackgroundHigh= 17e3;     % Altitude top to calcluate background    [m]
Def.BackgroundInd = 10;       % Pre-int bins to estimate background noise
Def.MaxTime       = 24*60*60; % Max time to run retrievals to         [sec]
Def.BinRange      = 2*37.5;   % Desired data range resolution         [m]
Def.MinRange      = 0;        % Start of retrievals                   [m]
Def.BlankRange    = 450;      % Altitude below which data is blanked  [m]
%% Defining Water Vapor retrieval options
% Standard retrieval
Op.WV             = Def;
Op.WV.Bootstrap   = false;    % Boolean to turn bootstrapping on
Op.WV.BootIters   = 30;       % Iterations to use when bootstraping
Op.WV.BinTime     = 4*60;     % Desired data time resolution          [sec]
Op.WV.GradFilt    = 1000;     % Count rate gradient to filter         [ ]
Op.WV.MaxRange    = 6.5e3;    % Max range to run retrievals to        [m]
Op.WV.MinRange    = 150;      % Start of retrievals                   [m]
Def.BlankRange    = 150;      % Altitude below which data is blanked  [m]
Op.WV.SmoothRange = 150;      % Desired smoothing range res           [m]
Op.WV.SmoothTime  = 4*60;     % Desired smoothing time res            [sec]
Op.WV             = MakeArrays(Op.WV);
% Low altitude retrieval
Op.WVLow             = Def;
Op.WVLow.BlankRange  = 90;      % Altitude below which data is blanked  [m]
Op.WVLow.Bootstrap   = false;   % Boolean to turn bootstrapping on
Op.WVLow.BootIters   = 30;      % Iterations to use when bootstraping
Op.WVLow.BinTime     = 4*60;    % Desired data time resolution        [sec]
Op.WVLow.BinRange    = 30;
Op.WVLow.GradFilt    = 1000;    % Count rate gradient to filter       [ ]
Op.WVLow.MaxRange    = 6.480e3; % Max range to run retrievals to      [m]
Op.WVLow.MinRange    = 0;       % Start of retrievals                 [m]
Op.WVLow.BlankCounts = 75;      % Altitude below which data is blanked[m]
Op.WVLow.SmoothRange = [];      % Desired smoothing range res         [m]
Op.WVLow.SmoothTime  = 6*60;    % Desired smoothing time res          [sec]
Op.WVLow             = MakeArrays(Op.WVLow);
%% Defining HSRL retrieval options
Op.HSRL             = Def;
Op.HSRL.Bootstrap   = false;  % Boolean to turn bootstrapping on
Op.HSRL.BootIters   = 30;     % Iterations to use when bootstraping
Op.HSRL.BinTime     = 2*60;   % Desired data time resolution          [sec]
Op.HSRL.BlankRange  = 250;    % Blank data from 0 to this range       [m]
Op.HSRL.MaxRange    = 8e3;    % Max range to run retrievals to        [m]
Op.HSRL.MinRange    = 150;    % Start of retrievals                   [m]
Op.HSRL.SmoothRange = 2*37.5; % Desired smoothing range res           [m]
Op.HSRL.SmoothTime  = 2*60;   % Desired smoothing time res            [sec]
Op.HSRL             = MakeArrays(Op.HSRL);
% Klett inversion parameters
Op.Klett.GuessK           = 1;
Op.Klett.GuessSigma       = 0.01;
Op.Klett.GuessLidarRatio  = 30;
%% Defining Temperature retrieval options
Op.Temp             = Def;
Op.Temp.Bootstrap   = false;  % Boolean to turn bootstrapping on
Op.Temp.BootIters   = 30;     % Iterations to use when bootstraping
Op.Temp.TempIter    = 50;     % Iterations for temperature conversion
Op.Temp.HSRLType    = 'On';   % Options: 'Py' (Python), 'On' (Calculated onboard),
                              % 'PyP' (Python Prefered if availible else onboard)
Op.Temp.Method      = 'LBL';  % Options: 'LBL' (line-by-line) or 'PCA'
Op.Temp.BinRange    = 4*37.5; % Desired data range resolution         [m]
Op.Temp.BinTime     = 5*60;   % Desired data time resolution          [sec]
Op.Temp.BlankBSC    = 5e-5;   % Blank data from 0 to this range       [m]
Op.Temp.MaxRange    = 6e3;    % Max range to run retrievals to        [m]
Op.Temp.MinRange    = 150;    % Start of retrievals                   [m]
Op.Temp.SmoothRange = 300;    % Desired smoothing range res           [m]
Op.Temp.SmoothTime  = 20*60;  % Desired smoothing time res            [sec]
Op.Temp             = MakeArrays(Op.Temp);
%% Overwritting settings if the user desires botstrapping
if BSOverwrite
    Op.WV.Bootstrap    = true;
    Op.WVLow.Bootstrap = true;
    Op.HSRL.Bootstrap  = true;
    Op.Temp.Bootstrap  = true;
    Op.Temp.HSRLType   = 'On';
end
%% Defining plotting options
Op.Plot.FontSize    = 18;     % Font size to be applied to all figures
Op.Plot.RB.MaxRange = 12e3;   % Max range to plot for field catalog   [m]
Op.Plot.RB.MinRange = 0;      % Min range to plot for field catalog   [m]
Op.Plot.RB.CAxis    = [1,6];  % Colorbar range for relative backscatter
Op.Plot.WV.MaxRange = 6e3;    % Max range to plot for field catalog   [m]
Op.Plot.WV.MinRange = 0;      % Min range to plot for field catalog   [m]
Op.Plot.WV.CAxis    = [0,10]; % Colorbar range for absolute humidity
%% Defining data to read 
Op.DataNames = {'Container';'Current';'Etalon';'HumiditySensor';'Laser';'MCS';
                'Power';'QuantumComposer';'Thermocouple';'UPS';'WeatherStation'}; 
end

function [S] = MakeArrays(S)
%
% Input: S: A strucutre containing user defined options
%
% Outputs: S: A stucutre modified with default time and range arrays
%
%% Adding arrays
S.MinTime   = S.BinTime./2;                      % Retrieval start  [sec]
S.Range     = S.MinRange:S.BinRange:S.MaxRange;  % Range array
S.TimeStamp = S.MinTime:S.BinTime:S.MaxTime;     % Time array
%% Forcing fields to be ordered
S = orderfields(S);
end
