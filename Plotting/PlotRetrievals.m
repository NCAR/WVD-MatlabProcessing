



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

%% Applying data masks and making my own 
Python = RecursivelyApplyMask(Python);

% Applying an optical depth filter to remove bad data above clouds
LidarRatio = 30;
Beta     = Python.BSCoefficient.Value;
Altitude = Python.BSCoefficient.Range;
% Calculating optical depth
Beta(isnan(Beta)) = 0;
OD =  double(2.*LidarRatio.*cumsum(Beta).*(Altitude(2)-Altitude(1)));
% Applying an optical depth filter
Python.BSCoefficient.Value(OD>0.25) = nan;


%% Plotting retrieved data
figure(FigNum);
set(gcf,'position',[1023,66,859,912],'color',[1,1,1])

%% Setting up Colormaps
Hot = flipud(colormap('hot'));

%% Plotting retrieved data
% Water vapor data
subplot(7,1,1:2);
pcolor(Python.Humidity.TimeStamp./60./60,Python.Humidity.Range./1e3,ApplyMask(Python.Humidity));
shading flat; colorbar; caxis([0,10]); colormap(gca,'parula')
xlim([0,24]);ylim([0,6]);
set(gca,'xtick',0:6:24,'ytick',0:1:6)
title([upper(erase(Options.System,'_')),' Retrievals (',Options.Date,')'])
AddPlotText(PO.TextXLoc,PO.TextYLoc,'Python Moisture [g/m^3]',PO.FontSize)
ylabel('Altitude [km]')

% Aerosol backscatter coefficient 
subplot(7,1,3:4);
pcolor(Python.BackRatio.TimeStamp./60./60,Python.BackRatio.Range./1e3,real(log10(ApplyMask(Python.BSCoefficient))));
shading flat; CB = colorbar; caxis([-8,-4]);
xlim([0,24]);ylim([0,6]);
set(gca,'xtick',0:6:24,'ytick',0:1:6)
AddPlotText(PO.TextXLoc,PO.TextYLoc,'Python Aerosol Backscatter Coefficient [m^{-1}sr^{-1}]',PO.FontSize)
ylabel('Altitude [km]')
colormap(gca,'jet')

% Temperature
subplot(7,1,5:6);
pcolor(Retrievals.Temperature.TimeStamp./60./60,Retrievals.Temperature.Range./1e3,Retrievals.Temperature.Smoothed-273.15);
shading flat; colorbar; caxis([-5,25])
xlim([0,24]);ylim([0,6]);
set(gca,'xtick',0:6:24,'ytick',0:1:6)
AddPlotText(PO.TextXLoc,PO.TextYLoc,'Temperature [^\circ C]',PO.FontSize)
ylabel('Altitude [km]')
colormap(gca,'parula')
DesiredPos = get(gca,'Position');

subplot(7,1,7)
plot(WS.TimeStamp,WS.Temperature);
xlim([0,24]);set(gca,'xtick',0:6:24)
AddPlotText(PO.TextXLoc,PO.TextYLoc,'Surface Temperature [^\circ C]',PO.FontSize)
xlabel('Hours UTC');ylabel('Temp [^\circ C]')
Current = get(gca,'Position');
set(gca,'Position',[Current(1:2),DesiredPos(3),Current(4)])

%% Format Figures
FormatFigures;

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

function [Data] = ApplyMask(Struct)
Data = Struct.Value;
Data(Struct.Mask==1) = nan;
end