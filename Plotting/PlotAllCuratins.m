% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created February 29, 2024

function [FigNum] = PlotAllCuratins(Data,Options,FigNum)
%
% Inputs: Data:    Structure containing all lidar data to be plotted
%         Options: Strucutre containing user defined plotting options
%         FigNum:  Desired figure number to plot (will be default if not
%                  used)
%
% Outputs: FigNum: Actual figure number plotted
%
%% Setting plot options
% Setting parameters for figure and panel size
PlotHeight   = 12;   % Total plot height
PlotWidth    = 16;   % Total plot width
SubplotCol   = 2;    % Number of subplot columns
SubplotRow   = 3;    % Number of subplot rows
LeftMargin   = 1;    % Margin on the left side of all subplots
RightMargin  = 1.75; % Margin on the right side of all subplots
TopMargin    = 1;    % Margin on the top side of all subplots
BottomMargin = 1.4;  % Margin on the bottom side of all subplots
SpaceRow     = 0.75; % Spacing between each row of subplots
SpaceCol     = 0.75; % Spacing between each column of subplots
FontSize     = 16;
% Determining where to place the figures
sub_pos=subplot_pos(PlotWidth,PlotHeight,LeftMargin,RightMargin,BottomMargin, ...
                    TopMargin,SubplotCol,SubplotRow,SpaceCol,SpaceRow);
% Flipping the figure list to fill from top to bottom
sub_pos = fliplr(sub_pos);
%% Checking inputs
if nargin == 2
    FigNum = FindCurrentFigure + 1;
end
%% Defining variables
Counter = 1;
%% Integrating and Background Subtracting Count Contours 
CountType = fields(Data);
for m=CountType'
    DataNeeded.(m{1}).TimeStamp = Data.(m{1}).TimeStamp.*60.*60;
    DataNeeded.(m{1}).Range = median(Data.(m{1}).RangeResolution,'omitnan').*1e-9.* ...
                          (1:1:size(Data.(m{1}).Data,2)).*3e8/2;
    DataNeeded.(m{1}).Counts = Data.(m{1}).Data';
end

% Plotting data with no additional binning
BinInfo.BinNum = [1;1];
BinInfo.BinDir = {'Range';'TimeStamp'};
% Defining the range and background settings
try
    Op.BinRange = DataNeeded.(CountType{1}).Range(2) - DataNeeded.(CountType{1}).Range(1);
catch
    Op.BinRange  = 37.5;
end
Op.BackgroundLow = 16e3;     % Altitude bottom to calcluate background [m]
Op.BackgroundHigh= 17e3;     % Altitude top to calcluate background    [m]
Op.MaxRange      = 12e3;
% Background subtracting 
Data2.BGSub = BGSubtractLidarData(DataNeeded,[],BinInfo,Op);

%% Plotting
figure(FigNum)
set(gcf,'position',[1933,150,1518,925],'color',[1,1,1])

for m={'WVOnline','WVOffline','O2OnlineMol','O2OfflineMol','O2OnlineComb','O2OfflineComb'}
    % Creating a grid on the time and range
    [x,y] = meshgrid(Data2.BGSub.(m{1}).TimeStamp./60./60, ...
                     Data2.BGSub.(m{1}).Range./1e3);
    % Plotting
    axes('position',sub_pos{Counter},'XGrid','off','XMinorGrid','off',...
                         'FontSize',FontSize,'Box','on','Layer','top');
    % Data2.BGSub.(m{1}).Counts(Data2.BGSub.(m{1}).Counts<0) = nan;
    try
        pcolor(x,y,Data2.BGSub.(m{1}).Counts)
    catch
        12
    end
    shading flat;
    if mod(Counter,2) == 1
        ylabel('Altitude [km]','FontSize',FontSize);
    end
    if Counter >= 5
        xlabel('UTC Hours','FontSize',FontSize)
    end
    if Counter == 6
        Position = get(gca,'position');
        Cb = colorbar;
        % Setting positions of last contour and colorbar
       set(gca,'position',Position)
       set(Cb,'position',[sub_pos{6}(1) + sub_pos{6}(3) + 0.02, sub_pos{6}(2),0.02, ...
                          sub_pos{2}(2) + sub_pos{2}(4) - sub_pos{6}(2)])
       ylabel(Cb,'Background Subtracted Counts','FontSize',FontSize)
    end
    xlim([0,24]); ylim([0,Op.MaxRange/1e3])
    set(gca,'xtick',0:4:24,'ytick',0:3:Op.MaxRange/1e3,... % Axis tick marks
        'FontSize',FontSize,'linewidth',2,'box','on',   ... % Axis appearance
        'xgrid','on','ygrid','on','layer','top',     ... % Adding gridlines     
        'ColorScale','log');                             % Adding gridlines
    clim([1e0,1e5])
    AddPlotText(0.02,0.9,m{1},FontSize)
    % Updating plot counter
    Counter = Counter + 1;
end
sgtitle(['MPD Background Subtracted Photon Counts (',Options.Date,')'],...
         'FontWeight','bold')

end

function AddPlotText(TextXLoc,TextYLoc,Text,FontSize)
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,Text,'Fontsize',FontSize,'Fontweight','bold')
end

function [ positions ] = subplot_pos(plotwidth,plotheight,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey)
 
    subxsize=(plotwidth-leftmargin-rightmargin-spacex*(nbx-1.0))/nbx;
    subysize=(plotheight-topmargin-bottommargin-spacey*(nby-1.0))/nby;
 
    for i=1:nbx
       for j=1:nby
 
           xfirst=leftmargin+(i-1.0)*(subxsize+spacex);
           yfirst=bottommargin+(j-1.0)*(subysize+spacey);
 
           positions{i,j}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth subysize/plotheight];
 
       end
    end
end