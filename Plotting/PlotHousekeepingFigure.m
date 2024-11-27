

function [FigNum] = PlotHousekeepingFigure(Data,Options,FigNum)
%
%
%
%
%% Checking inputs
if nargin == 2
    FigNum = FindCurrentFigure + 1;
end
%%
figure(FigNum);
set(gcf,'position',[1023,66,859,912],'color',[1,1,1])
%% Plotting Options
PO.Hue      = {[0 0 0];[1 0 0];[0 0 1];[.8 0 .8];[0 1 1];[0 .5 0];};           % k,r,b,c,m,g
PO.TextXLoc = 0.01;
PO.TextYLoc = 0.9;
PO.FontSize = 12;
PO.XLims    = [0,24]; 
PO.XTick    = 0:6:24;
PO.MaxSubs  = 8*4;
Subs        = [];
%% Temperatures
Subs = GeneralPlotter(Data,{'UPS';'WeatherStation';{'Thermocouple','OpticalBench'};{'Thermocouple','HVACSource'};'HumiditySensor';{'Thermocouple','WVEtalonHeatSink'}},...
                      {'Temperature';'Temperature';'Temperature';'Temperature';'IntTemp';'Temperature'},...
                      {{'Temperature';'[^oC]'}},{1:4},{[-45,45]},{[]},Subs,PO); 

% Add second transparent axis on top of the last with pressure
Subs = CopyPlotter(Data,{{'Thermocouple','HSRLOven'}},{'Temperature'},...
                      {{'HSRL Cell';'[^oC]'}},{[0,120]},{[]},Subs,PO);

title([upper(erase(Options.System,'_')),' Housekeeping Parameters (',Options.Date,')'])
% Add second transparent axis on top of the last with HSRL cell temperature

%% Environmental Data
Subs = GeneralPlotter(Data,{'WeatherStation';'HumiditySensor'},...
                      {'RelativeHumidity';'RelativeHumidity'},...
                      {{'Relative';'Hum. [%]'}},{5:8},{[0,100]},{[]},Subs,PO);
% Add second transparent axis on top of the last with pressure
Subs = CopyPlotter(Data,{'WeatherStation'},{'Pressure'},...
                      {{'Pressure';'[mbar]'}},{[0,1100]},{[]},Subs,PO);
%% Laser Wavelength/Current/Seed Power
Subs = LaserPlotter(Data,{'TimeSeries','Laser'},{9:12;13:16;17:20},...
                    {'WaveDiff';'Current';'SeedPower'},[1000,1000,1],...
                    {{'Laser \lambda';'Deviation [pm]'},{'Laser Current';'[mA]'},{'Laser Seed';'Power [dBm]'}},...
                    {[-10,10],[0,200],[-25,5]},{[-.2,.2],5,1},Subs,PO);
%% Laser Amplified Power
Subs = LaserPlotter(Data,{'TimeSeries','Power'},{21:24},{'LaserPower'},1e-5,...
                       {{'Laser Power';'x 10^5 [arb]'}},{[]},{[]},Subs,PO);
%% Etalon Temperature Deviation
Subs = LaserPlotter(Data,{'TimeSeries','Etalon'},{25:28},{'TempDiff'},1,...
                      {{'Etalon Temp';'Deviation [^oC]'}},{[-5,5]},...
                      {[-0.05,.05]},Subs,PO);
%% Current 
Subs = GeneralPlotter(Data,{{'Current','SystemInput'};{'Current','HVAC'};{'Current','WallPowerStrip'};{'Current','WallHeaters'};{'Current','WindowHeaterFan'}},...
                      {'Current';'Current';'Current';'Current';'Current'},...
                      {{'Current';'[A]'}},{29:32},{[0,30]},{[0,15]},Subs,PO);
%% Moving the y-axes from side to side for readability 
for m=3:1:size(Subs,2)
    if mod(m,2) == 0
        set(Subs(m),'YAxisLocation','right')
    else
        set(Subs(m),'YAxisLocation','left')
    end
end
end

function [Subs] = CopyPlotter(Data,Types2Plot,Var2Plot,SubLabel,YMaxLim,YMinLim,Subs,PO)
%
% Inputs: 
%
%
%
%
%
%
%
%% Plotting constants
Color2Use = [75,142,173]./255;
%% Recreating the details of the axis to copy
CurrentAxis = gca;
Position = get(CurrentAxis,'position');
axes('Position',Position,'YAxisLocation','right','Color','none')
set(gca,'ycolor',Color2Use)
%% Looking for availible data and plotting if there
for m=1:1:size(Types2Plot,1)
   [IsField,TempData] = RecursivelyCheckIsField(Data,cellflat({'TimeSeries',Types2Plot{m,1}}));
   if IsField
       hold on; plot(TempData.TimeStamp,TempData.(Var2Plot{m}),'Color',Color2Use)
   end
end
%% Formatting subplot
Subs = FormatGeneral('',[],YMaxLim,YMinLim,SubLabel,Subs,PO);
% Putting original axis on top
uistack(CurrentAxis,'top')
end

function [Subs] = GeneralPlotter(Data,Types2Plot,Var2Plot,SubLabel,Subplots,YMaxLim,YMinLim,Subs,PO)
%
% Inputs: 
%
%
%
%
%
%
%
%%
Label    = '';
%% Looking for availible data and plotting if there
n=1;
for m=1:1:size(Types2Plot,1)
   [IsField,TempData] = RecursivelyCheckIsField(Data,cellflat({'TimeSeries',Types2Plot{m,1}})); 
   if iscell(Types2Plot{m,1})
       Label = [Label,'\color[rgb]{',num2str(PO.Hue{m}),'} ',Types2Plot{m}{end},'\color{black}, ']; %#ok<AGROW>
   else
       Label = [Label,'\color[rgb]{',num2str(PO.Hue{m}),'} ',Types2Plot{m},'\color{black}, ']; %#ok<AGROW>
   end
   if IsField
       subplot(PO.MaxSubs,1,Subplots{n}); hold on;
       plot(TempData.TimeStamp,TempData.(Var2Plot{m}),'Color',PO.Hue{m})
   end
end
%% Formatting subplot
Subs = FormatGeneral(Label,Subplots,YMaxLim,YMinLim,SubLabel,Subs,PO);
end

function [Subs] = LaserPlotter(Data,Data2Find,Subplots,SubVars,Multiple,SubLabel,YMaxLim,YMinLim,Subs,PO)
%
% Inputs: Data:      Data strucutre containing all of the processed data
%         Data2Find: Substructure to look for in Data. 
%         Subplots:  A cell of elements defining which subplots to use when 
%                    creating a subplot
%         SubVars:   A cell of variables to plot coming from the subplot
%                    pulled from the Data with Data2Find
%         Multiple:  An array of constant offsets to use to scale the
%                    element to be plotted
%         SubLabel:  A call of strings used to put on each subplot
%
%%
Label    = '';
%% Looking for availible data and plotting if there
[IsField,TempData] = RecursivelyCheckIsField(Data,Data2Find);
% If the data is availible in the general data structure, plot it
if IsField
    % Determining what elements are availible and removing unknown data
    Av = fieldnames(TempData);
    Av(ismember(Av,'Unknown')) = []; 
    % Plotting all data
    for n=1:1:size(Subplots)
        Label    = '';
        for m=1:1:size(Av,1)
            % Only plot Amplifier data if it is current
            if strcmp(Av{m}(end-2:end),'Amp') == 0 || strcmp(SubVars{n},'Current')
                subplot(PO.MaxSubs,1,Subplots{n}); hold on;
                plot(TempData.(Av{m}).TimeStamp,TempData.(Av{m}).(SubVars{n}).*Multiple(n),'Color',PO.Hue{m})
                Label = [Label, '\color[rgb]{',num2str(PO.Hue{m}),'}',Av{m},'\color{black}, ']; %#ok<AGROW>
            end
        end
        L{n} = Label;
    end
end
%% Formatting subplot
Subs = FormatGeneral(L,Subplots,YMaxLim,YMinLim,SubLabel,Subs,PO);
end

function [Subs] = FormatGeneral(Label,Subplots,YMaxLim,YMinLim,SubLabel,Subs,PO)
%
%
%
%
%%
ExtraPlot = false; UseGrid = 'on';
if isempty(Subplots)
    ExtraPlot = true;
    Subplots = {-1};
    UseGrid     = 'off';
end
%%
for n=1:1:length(Subplots)
    if ~ExtraPlot
        Subs(end+1) = subplot(PO.MaxSubs,1,Subplots{n});  %#ok<AGROW>
    end
    hold on;
    % Checking to see if the current y-axis limits are reasonable
    if isempty(YMaxLim{n}) ~= 1; ControlYMaxLims(YMaxLim{n}); end
    if isempty(YMinLim{n}) ~= 1
        if length(YMinLim{n}) == 1
            YLimits = ylim;
            ControlYMinLims([YLimits(1)-YMinLim{n},YLimits(2)+YMinLim{n}]);
        else
            ControlYMinLims(YMinLim{n});
        end
    end
    % Defining the base format of the subplot
    set(gca,'Fontsize',PO.FontSize,'Xlim',PO.XLims,'xtick',PO.XTick,...
        'box','on','xgrid',UseGrid,'ygrid',UseGrid,'Color','none',...
        'LineWidth',1.5)
    % Labeling the y-axis
    ylabel(SubLabel{n},'Fontsize',PO.FontSize,'FontWeight','bold')
    % Removing the x-axis labels if not the bottom subplot
    if max(Subplots{n}) ~= PO.MaxSubs; set(gca,'xticklabel',{}); end
    % Adding text to define subplot contents (legend)
    if iscell(Label)
        L = Label{n}(1:end-2);
    else
        L = Label(1:end-2);
    end
    AddPlotText(PO.TextXLoc,PO.TextYLoc,L,PO.FontSize)
end
end

function ControlYMaxLims(MaxLimits)
%
% Inputs: A 2 element array used to specify the user maximum y-axis limits
%
%%
YLimits = ylim;
if YLimits(1) < MaxLimits(1); YLimits(1) = MaxLimits(1); end
if YLimits(2) > MaxLimits(2); YLimits(2) = MaxLimits(2); end
if YLimits(1) > YLimits(2); YLimits(1) = YLimits(2) - 10; end
ylim(YLimits);
end

function ControlYMinLims(MinLimits)
%
% Inputs: A 2 element array used to specify the user minimum y-axis limits
%
%%
YLimits = ylim;
if YLimits(1) > MinLimits(1); YLimits(1) = MinLimits(1); end
if YLimits(2) < MinLimits(2); YLimits(2) = MinLimits(2); end
ylim(YLimits);
end

function AddPlotText(TextXLoc,TextYLoc,Text,FontSize)
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,Text,'Fontsize',FontSize,'Fontweight','bold')
end

function out = cellflat(celllist)
% CELLFLAT is a helper function to flatten nested cell arrays. 
% 
% CELLFLAT(celllist) searches every cell element in cellist and put them on
% the top most level. Therefore, CELLFLAT linearizes a cell array tree
% structure. 
%
% Example: cellflat({[1 2 3], [4 5 6],{[7 8 9 10],[11 12 13 14 15]},{'abc',{'defg','hijk'},'lmnop'}}) 
% 
% Output: 
%Columns 1 through 7
%     [1x3 double]    [1x3 double]    [1x4 double]    [1x5 double]    'abc'    'defg'    'hijk'
%   Column 8 
%     'lmnop'
%
% cellflat(({{1 {2 3}} 'z' {'y' 'x' 'w'} {4 @iscell 5} 6}) )
% Output: 
% [1]    [2]    [3]    'z'    'y'    'x'    'w'    [4]    @iscell    [5]    [6]
%
% Version: 1.0
% Author: Yung-Yeh Chang, Ph.D. (yungyeh@hotmail.com)
% Date: 12/31/2014
% Copyright 2015, Yung-Yeh Chang, Ph.D.
% See Also: cell
if ~iscell(celllist)
    error('CELLFLAT:ImproperInputAugument','Input argument must be a cell array');
end
out = {};
for idx_c = 1:numel(celllist)
    if iscell(celllist{idx_c})
        out = [out cellflat(celllist{idx_c})];
    else
        out = [out celllist(idx_c)];
    end
end
end
