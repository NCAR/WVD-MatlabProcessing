


function [Labels,FigNum] = PlotStatusFigure(Data,RawData,Options,FigNum)
%
%
%
%
%
%% Plotting Constants
Bounds.QError       = [20,100];     % Elements in queues before (warning, error)
Bounds.LLError      = [1e-4,5e-4];  % Stability of lasers before warnings         [nm]
Bounds.LLSeedStable = [0.5,1];      % Stability of laser seeds before warnings    [dBM]
Bounds.LLSeedLow    = [20,25];      % Minimum value of seed laser power times -1  [dBM]
Bounds.EtalonStable = [0.05,0.25];  % Stability of etalons before warnings        [C]
%% Checking inputs
if nargin == 3
    FigNum = FindCurrentFigure + 1;
end
%% Setting label information
LabelInfo = DefineLabelTypes(Options.System);
%% Making the labels
[Labels,Last] = MakeLabels(LabelInfo.Children(:,1),{'Responding';'Queues';'Data'});   % Children
A = cellfun(@(X) strcat(X,{'Locked';'SeedPower'}),LabelInfo.Lasers,'Uni',false);
[Labels,Last] = MakeLabels({''},cat(1,A{:}),Labels,Last);                             % Lasers
[Labels,Last] = MakeLabels({''},strcat(LabelInfo.Etalons,{'Locked'}),Labels,Last);    % Etalons 
[Labels,Last] = MakeLabels({'UPS'},LabelInfo.UPS,Labels,Last);                        % UPS
[Labels, ~  ] = MakeLabels({''},strcat(LabelInfo.Environmental,{'Temp'}),Labels,Last);% Environmental 
clear A Last
%% Pre-allocating data
Time         = Options.TimeGrid1d;
Contour      = nan(size(Time,1),size(Labels,1));
DefaultData  = -1*ones(size(Time));
GLTO         = -1;     % Global Last Time Observed
%% Unpacking the MPD data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Children %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:1:size(LabelInfo.Children,1)
   % Default data if if can not be found
   DataAvail = DefaultData; Responding = DefaultData; Queues = DefaultData;
   % Checking for end of data (children with 2 data types, pick reliable one)
   if strcmp(LabelInfo.Children{m,2},'WavelengthLocking'); CheckData = 'Laser';
   elseif strcmp(LabelInfo.Children{m,2},'MCS');           CheckData = 'Power';
   else;                                                   CheckData = LabelInfo.Children{m,2};
   end
   % Checking if data is availible
   [~,TempData] = RecursivelyCheckIsField(RawData,CheckData);
   if isempty(TempData) 
       DataAvail(:) = 2;
   else
       [Raw,GLTO] = DetermineDataAvailibility(TempData,Time,Options,GLTO,true);
       if isempty(Raw) == 0
           % Determining where time stamps are bad
           DataAvail(~isnan(Raw)) = 0;
           DataAvail(isnan(Raw)) = 2;
       end
   end
   % Checking if data exists
   [IsField,TempData] = RecursivelyCheckIsField(Data, {'TimeSeries','Container',LabelInfo.Children{m,2}});
   if IsField
       % Checking if the child is responding
       if isfield(TempData,'Status')
           Responding(TempData.Status==1) = 0;
           Responding(TempData.Status==0 | isnan(TempData.Status)) = 2;
       end
       % Check if the child queues are overflowing
       if isfield(TempData,'CQueueEl') && isfield(TempData,'RQueueEl')
           Queues = max([CheckDataStatus(Queues,TempData.RQueueEl,Bounds.QError),...
                         CheckDataStatus(Queues,TempData.CQueueEl,Bounds.QError)],[],2);
           Queues(isnan(TempData.RQueueEl)|isnan(TempData.CQueueEl)) = 2;
           % Adding check: if not responding...queue data is invalid
           Queues(Responding == 2) = 2;
       end
   end
   % Filling the data contour
   Contour(:,3.*(m-1)+(1:3)) = [Responding,Queues,DataAvail];
end
LastContour = 3.*m;
clear CheckData DataAvail Responding Queues IsField TempData m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lasers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:1:size(LabelInfo.Lasers,1)
    % Default data if if can not be found
    Locked = DefaultData; Stable = DefaultData; Low = DefaultData; 
    % Checking if data exists
    [IsField,TempData] = RecursivelyCheckIsField(Data, {'TimeSeries','Laser',LabelInfo.Lasers{m}});
    % If data exists, determine what its color coding should be
    if IsField
        % Checking for the end of the data
        [~,GLTO] = DetermineDataAvailibility(TempData,Time,Options,GLTO,false);
        % Checking if the lasers are locked
        Locked = CheckDataStatus(Locked,TempData.WaveDiff,Bounds.LLError);
        % Checking if the seeds are stable 
        if sum(~isnan(TempData.SeedPower)) == 0
            FirstGood = -70;
        else
            FirstGood = TempData.SeedPower(~isnan(TempData.SeedPower));
        end
        Stable = CheckDataStatus(Stable,TempData.SeedPower - FirstGood(1),Bounds.LLSeedStable);
        Low    = CheckDataStatus(Low,TempData.SeedPower,Bounds.LLSeedLow);
        Stable(Low == 1) = 1; Stable(Low == 2) = 2;
    end
    % Filling the data contour
    Contour(:,LastContour + (m-1)*2 + (1:2)) = [Locked,Stable];
end
LastContour = LastContour + 2*m;
clear Locked Stable IsField TempData m Raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Etalons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:1:size(LabelInfo.Etalons,1)
    % Default data if if can not be found
    Locked = DefaultData;
    % Checking if data exists
    [IsField,TempData] = RecursivelyCheckIsField(Data, {'TimeSeries','Etalon',LabelInfo.Etalons{m}});
    % If data exists, determine what its color coding should be
    if IsField
        % Checking for the end of the data
        [~,GLTO] = DetermineDataAvailibility(TempData,Time,Options,GLTO,false);
        Locked = CheckDataStatus(Locked,TempData.TempDiff,Bounds.EtalonStable);
    end
    % Filling the data contour
    Contour(:,LastContour + m) = Locked;   
end
LastContour = LastContour + m;
clear Locked IsField TempData m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:1:1
    Battery    = DefaultData;
    OnBattery  = DefaultData;
    BatteryLow = DefaultData;
    % Checking if data exists
    [IsField,TempData] = RecursivelyCheckIsField(Data, {'TimeSeries','UPS'});
    % If data exists, determine what its color coding should be
    if IsField
        % Checking for the end of the data
        [~,GLTO] = DetermineDataAvailibility(TempData,Time,Options,GLTO,false);
        % If data exists, determine what its color coding should be
        if isfield(TempData,'BatteryInUse')
            OnBattery(TempData.BatteryInUse == 1) = 0;
            OnBattery(TempData.BatteryInUse == 0) = 2;
        end
        if isfield(TempData,'BatteryLow')
            BatteryLow(TempData.BatteryLow == 0) = 0;
            BatteryLow(TempData.BatteryLow == 1) = 2;
        end
        if isfield(TempData,'BatteryReplace') && isfield(TempData,'BatteryNominal')
            Battery(TempData.BatteryReplace == 0 & TempData.BatteryNominal == 1) = 0;
            Battery(TempData.BatteryReplace == 1 | TempData.BatteryNominal == 0) = 1;
        end
    end
    % Filling the data contour
    Contour(:,LastContour + (1:3)) = [Battery,OnBattery,BatteryLow];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thermocouples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Removing data not yet taken
% I should check here to see if the global last time observed is low
% because there is no data at the end of the day or because it could not
% have possibly been taken yet
if GLTO > 0 && (now-1 < datenum(Options.Date,'yyyymmdd'))
    % Final check...the absolute Global Last Observed Time should be right
    % now. If that is less than the GLTO, there was a file issue
    TimeOfDay = (now-datenum('20210715','yyyymmdd'))*24;
    if GLTO > TimeOfDay && TimeOfDay>0
        GLTO = TimeOfDay;
    end
    Contour((Time > GLTO),:) = nan;
end
%% Plotting all data
% Setting up the elements needed to plot and label the axis
StatusElements = size(Labels,1);
[X,Y] = meshgrid(Time,0.5:StatusElements-0.5);
% Plotting the contour on the left axis
figure(FigNum);
yyaxis left
StandardAxisPrep(StatusElements,Labels(:,2))
% Plotting the names on the right axis
yyaxis right
pcolor(X,Y,Contour');
StandardAxisPrep(StatusElements,Labels(:,3));
% Labeling the x-axis and titles
xlabel('Time [UTC]'); 
title({[upper(erase(Options.System,'_')),' Status (',Options.Date,')'];
    ['\color{black}(\color[rgb]{0.7,0.7,0.7}data missing\color{black},',...
     ' \color[rgb]{0,0.5,0}nominal\color{black}, \color[rgb]{0.8,0.8,0}warning',...
     '\color{black}, \color[rgb]{0.75,0.00,0.00}fault\color{black})']})
%% Formatting the figure
FormatStatusFigure(gca)
end

function [LabelInfo] = DefineLabelTypes(Type)
%
%
%
%
%% Labels for the running labview programs ({matlab label,labview label})
LabelInfo.Children       = {'MCS'    ,'MCS';           
                            'LL'     ,'WavelengthLocking';
                            'WS'     ,'WeatherStation';
                            'HK'     ,'Thermocouple';
                            %'Cur'    ,'CurrentSensor';
                            'Hum'    ,'HumiditySensor';
                            'UPS'    ,'UPS';
                            'Clock'  ,'QuantumComposer'};
%% Settings for installed hardware
switch Type
%     case 'mpd_02'
%         LabelInfo.Lasers        = {'WVOnline';'WVOffline';'HSRL'};
%         LabelInfo.Etalons       = {'WVEtalon';'HSRLEtalon';'HSRLEtalon2'};
%         LabelInfo.Environmental = {'WVEtalonHeatSink';'RbEtalonHeatSink';'Bench';'RbOven';'Humidity'};
    case 'mpd_05'
        LabelInfo.Lasers        = {'WVOnline';'WVOffline';'O2Online';'O2Offline'};
        LabelInfo.Etalons       = {'WVEtalon';'O2Etalon'};
        LabelInfo.Environmental = {'WVEtalonHeatSink';'O2EtalonHeatSink';'Bench';'KOven';'Humidity'};
    otherwise
        LabelInfo.Lasers        = {'WVOnline';'WVOffline'};
        LabelInfo.Etalons       = {'WVEtalon'};
        LabelInfo.Environmental = {'WVEtalonHeatSink';'Bench';'Humidity'};
end
LabelInfo.UPS            = {'BatteryStatus';'OnBattery';'BatteryLow'};

end

function [Data] = CheckDataStatus(Data,Evaluator,Bounds)
%
%
%
%% Checking how data compares to [Good/Warning/Bad] thresholds
Data(abs(Evaluator) <  Bounds(2)) = 1;
Data(abs(Evaluator) <  Bounds(1)) = 0;
Data(abs(Evaluator) >  Bounds(2)) = 2;
end

function [Labels,Last] = MakeLabels(FirstItem,Status,Labels,Last)
%
%
%
%
%
%
%
%%
if nargin == 2
    Labels = {};
    Last   = [];
end
%% Creating the labels and determining what side they should be on
for m=1:1:size(FirstItem,1)
    for n=1:1:size(Status,1)
        Labels{end+1,1} = [FirstItem{m},' ',Status{n}];  %#ok<*AGROW>
        Side = DetermineSide(m,Last);
        if strcmp(Side,'Left')
            Labels{end,2}   = Labels{end,1};
            Labels{end,3}   = '';
        else
            Labels{end,2}   = '';
            Labels{end,3}   = Labels{end,1};
        end
    end
end
Last = Side;
end

function [Return] = DetermineSide(m,Last)
%
%
%
%
%
%
%%
CheckVar = 0;
if isempty(Last) || strcmp(Last,'Right')
    CheckVar = 1;
end
if mod(m,2) == CheckVar
    Return = 'Left';
else
    Return = 'Right';
end
end

function [Raw,GLTO] = DetermineDataAvailibility(TempData,Time,Options,GLTO,CalcDA)
%
% InputsL Temp:    
%         Time:    
%         Options: 
%         GLTO:    Global last time observed
%         CalcDA:  Bool to check whether to calculate data availibility 
%                  
% Outputs: Raw:    
%          GLTO:   
%
%% Looking for data breaks
try TempData = rmfield(TempData,'Type'); catch; end
TempData = RecursivelyIdentifyBreaks(TempData,Options.BreakSize);
%%
[~,Raw] = RecursivelyCheckIsField(TempData,'TimeStamp');
if isempty(Raw) == 0
    % Forcing time stamps to be monotonic (because they come from raw data)
    Raw = ForceMonotonicTimeStamps(Raw);
    % Finding data breaks if they exist to apply later
    FieldNames = fieldnames(rmfield(TempData,'TimeStamp'));
    DataCheck = TempData.(FieldNames{1,1})(:,1);
    % Check if first time stamp is causing a full day shift
    for m=1:1:5  % Only moving back 5 days (more than 1 would be weird)
        if mean(Raw,'omitnan') > 24 % Shifting back a full day to check
            Raw = Raw - 24;
        else 
            break
        end
    end
    % Determining the max data time
    GLTO = max([max(Raw),GLTO]);
    % Interpolating raw time stamps to known time grid if we are looking at
    % data availibility. Otherwise this step is prone to error and unneeded
    if CalcDA
        DataCheck = interp1(Raw,DataCheck,Time);
        Raw = Time;
        Raw(isnan(DataCheck)) = nan;
    end
end
end

function StandardAxisPrep(StatusElements,Labels2Use)
%
%
%
%
%
%% 
MinX     = 0;
MaxX     = 24;
FontSize = 12;
%% Setting axis parameters
xlim([MinX,MaxX]); 
ylim([0.5,StatusElements+0.5]); 
set(gca,'Ydir','reverse');
set(gca,'xtick',MinX:3:MaxX,'ytick',1:StatusElements,'yticklabel',Labels2Use);
set(gca,'ycolor','k','fontsize',FontSize)
shading flat;
%% Setting up the colormap
Colormap = [0.70,0.70,0.70;       % White  = -1
            0.00,0.50,0.00;       % Green  = 0
            0.80,0.80,0.00;       % Yellow = 1
            0.75,0.00,0.00];      % Red    = 2
colormap(Colormap);
caxis([-1.5,2.5])
end

function FormatStatusFigure(CurrentFigure)
%
% Inputs:
%
% Outputs: none
%
%% Figure style
FontSize        = 14;
ILW             = 0.25;                % Inner Line Width
OLW             = 4;                   % Outter Line Width
FigureHeight    = 900;
FigureWidth     = 900;
MenuHeight      = 90;
ScreenSize      = get(0,'ScreenSize');
%% Setting figure parameters
% Sizing and coloring the figure
set(CurrentFigure,'Color',[1 1 1],'FontSize',FontSize); hold on;
set(gcf,'Color',[1 1 1],'Position',[10 (ScreenSize(4)-FigureHeight-MenuHeight)  FigureWidth FigureHeight]);
% Bringing in the plot to make y-labels less squished
Position        = get(gca,'position');
HorShift        = 0.03;
VerShift        = 0.01;
Position(1)     = Position(1) + HorShift;
Position(2)     = Position(2) - VerShift;
Position(3)     = Position(3) - 2.*HorShift;
Position(4)     = Position(4) + 2.*VerShift;
set(gca,'position',Position)
%% Line Locations
XLims = xlim;
YLims = ylim;
XTick = get(gca,'xtick');
YTick = YLims(1):1:YLims(2);
Pts   = 1000;
%% Plotting lines
% Vertical
for m=1:1:size(XTick,2)
    TopPlot(XTick(m)*ones(1,Pts),linspace(YLims(1),YLims(2),Pts),'k',':',ILW)
end
% Horizontal 
for m=1:1:size(YTick,2)-1
    TopPlot(linspace(XLims(1),XLims(2),Pts),YTick(m)*ones(1,Pts),'k','--',ILW)
end
% Outside box 
for m=1:1:2    
    TopPlot(linspace(XLims(1),XLims(2),Pts),YLims(m)*ones(1,Pts),'k','-',OLW)
    TopPlot(XLims(m)*ones(1,Pts),linspace(YLims(1),YLims(2),Pts),'k','-',OLW)
end
end
