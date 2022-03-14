% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created March, 2022

function FigNum = PlotWVQuicklook(WV,Options,Op)
%
% Inputs: WV:
%         Options: Options structure for just WV
%         Op: Full options structure
%
% Outputs: none
%
%% Setting up needed variables
Sc   = get(0,'ScreenSize');
Date = datestr(datenum(Op.Date,'yyyymmdd'), 'dd mmm yyyy');
%% Loading colormap 
RBColormap = importdata('NCAR_C_Map.mat');
%% Formatting the figure as a whole
FigNum = FindCurrentFigure + 1;
figure(FigNum)
set(gcf,'Position',[Sc(3)/20 Sc(4)/10 Sc(3)/1.5 Sc(4)/1.5]);
%% Plotting Relative Backscatter
subplot(2,1,1)
pcolor(WV.RB.TimeStamp./60./60,WV.RB.Range./1e3,real(log10(WV.RB.Value*5)));
colormap(gca,RBColormap)
% Formatting the subplot
cb = FormatAxis(Date, Options.RB, 'Relative Backscatter (C/ns km^2)');
% Reformatting colorbar ticks to make user know it is logrithmic
A = get(cb,'yticklabel');
C = get(cb,'ytick');
for m=1:1:size(A,1); B{m,1} = ['10^{',A{m},'}']; end      %#ok<AGROW>
set(cb,'ytick',C,'yticklabel',B)
%% Plotting Water Vapor
% Applying the water vapor data mask
WV.Smoothed2(WV.Mask) = nan;
% Plotting
subplot(2,1,2)
pcolor(WV.TimeStamp./60./60,WV.Range./1e3,WV.Smoothed2);
colormap(gca,'jet')
% Formatting the subplot
FormatAxis(Date, Options.WV, 'Water Vapor (g/m^3)');
%% Formatting figure
FormatFigures;
end

function CB = FormatAxis(Date, Op, Name)
%
%
%
%
%
%% Formatting the pcolor type
shading interp
%% Labeling plot
xlabel('Time (UTC)');
ylabel('Height (km, AGL)');
title({[Date,' ',Name]});
%% Formatting axes
set(gca,'xtick',0:2:24)
%% Setting axis bounds
xlim([0,24]); ylim([Op.MinRange,Op.MaxRange]./1e3);
%% Formatting the colorbar
caxis(Op.CAxis); CB = colorbar('EastOutside');
end

