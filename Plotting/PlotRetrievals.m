% Written By: Robert Stillwell
% Written For: NCAR

function [FigNum] = PlotRetrievals(Retrievals,Python,Options,WS,FigNum)
%
%
%
%
%
%% Checking inputs
if nargin ~= 5
    FigNum = FindCurrentFigure + 1;
end
%% Plotting Options
PO.TextXLoc = 0.01;
PO.TextYLoc = 0.9;
PO.FontSize = 16;

%% Applying optical depth filter to HSRL data
Retrievals.HSRL.Smoothed(Retrievals.HSRL.OD>3.0) = nan;
Retrievals.HSRL.Smoothed(Retrievals.HSRL.Mask) = nan;

%% Plotting retrieved data
figure(FigNum);
set(gcf,'position',[1023,66,859,912],'color',[1,1,1])
%% Plotting retrieved data
% Water vapor data
subplot(7,1,1:2);
pcolor(Retrievals.WaterVapor.TimeStamp./60./60,Retrievals.WaterVapor.Range./1e3,Retrievals.WaterVapor.Smoothed2)
FormatAxis([0,10],CM_viridis(64),Options,PO,'Absolute Humidity [g/m^3]','Retrievals',[],[0,24],[0,6]);
% Aerosol backscatter coefficient 
subplot(7,1,3:4);
pcolor(Retrievals.HSRL.TimeStamp./60./60,Retrievals.HSRL.Range./1e3,real(log10(Retrievals.HSRL.Smoothed)))
CB = FormatAxis([0,2],flipud(CM_magma(64)),Options,PO,'Backscatter Ratio [unitless]',[],[],[0,24],[0,6]);
% Temperature
subplot(7,1,5:6);
pcolor(Retrievals.Temperature.TimeStamp./60./60,Retrievals.Temperature.Range./1e3,Retrievals.Temperature.Smoothed-273.15);
FormatAxis([-5,25],CM_cividis(64),Options,PO,'Temperature [^\circ C]',[],[],[0,24],[0,6]);
DesiredPos = get(gca,'Position');
% Surface station
subplot(7,1,7)
plot(WS.TimeStamp,WS.Temperature);
FormatAxis([],[],Options,PO,'Surface Temperature [^\circ C]',[],'Temp [^\circ C]',[0,24],ylim);
xlabel('Hours UTC');
Current = get(gca,'Position');
set(gca,'Position',[Current(1:2),DesiredPos(3),Current(4)])

%% Format Figures
FormatFigures(gcf);

%% Formatting colorbars
ColorbarTicks = get(CB,'yticklabel');
ColorbarTicksNew = ColorbarTicks;
for a=1:1:size(ColorbarTicks,1)
    ColorbarTicksNew{a} = ['10^{',ColorbarTicks{a},'}'];
end
set(CB,'yticklabel',ColorbarTicksNew)

end

function AddPlotText(TextXLoc,TextYLoc,Text,FontSize)
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,Text,'Fontsize',FontSize,'Fontweight','bold')
end

function [CB] = FormatAxis(CBounds,CMap,Op,PO,Text,Title,YLabel,XLim,YLim)
% Setting up the color and colorbar
if ~isempty(CBounds)
    shading flat; CB = colorbar; caxis(CBounds); colormap(gca,CMap)
end
% Setting up the x/y axis and bounds
xlim(XLim);ylim(YLim);
set(gca,'xtick',XLim(1):6:XLim(2))
% Labeling
if~isempty(Title)
    title([upper(erase(Op.System,'_')),' ',Title,' (',Op.Date,')'])
end
AddPlotText(PO.TextXLoc,PO.TextYLoc,Text,PO.FontSize)
if isempty(YLabel)
    YLabel = 'Altitude [km]';
end
ylabel(YLabel)
end
