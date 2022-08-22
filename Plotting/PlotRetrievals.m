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
% Approximating backscatter coefficient
LidarRatio = 30;
[BetaM,~] = RayleighBackscatterCoeff(770.1085e-9,Retrievals.HSRL.PGuess.*1013.25,Retrievals.HSRL.TGuess);
BSR = Retrievals.HSRL.SmoothedV2;
Beta = BetaM.*(BSR-1);
Beta(isnan(Beta)) = 0;
% Calculating optical depth from backscatter coefficient
OD = double(2.*LidarRatio.*cumsum(Beta).*(Retrievals.HSRL.Range(2)-Retrievals.HSRL.Range(1)));
% Applying an optical depth filter
Retrievals.HSRL.SmoothedV2(OD>3.0) = nan;

%% Copying HSRL and WV masks for Temperature
TempMask = isnan(Python.MPD.BSCoefficient.Value) | isnan(Python.MPD.Humidity.Mask);
% Interpolating Python grid to temperature grid
[X,Y] = meshgrid(Python.MPD.BackRatio.TimeStamp, Python.MPD.BackRatio.Range);
[x,y] = meshgrid(Retrievals.Temperature.TimeStamp,Retrievals.Temperature.Range);
TempMaskInt = ceil(interp2(X,Y,double(TempMask),x,y));
% Applying mask
Retrievals.Temperature.Smoothed(TempMaskInt==1) = nan;
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
pcolor(Retrievals.HSRL.TimeStamp./60./60,Retrievals.HSRL.Range./1e3,real(log10(Retrievals.HSRL.SmoothedV2)))
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

function [Beta,BetaTotal] = RayleighBackscatterCoeff (Lambda,Press,Temp)
%
% Inputs: Lambda:     Laser wavelength                   [meters]
%         Press:      Atmospheric pressure               [millibar]
%         Temp:       Atmospheric temperature            [Kelvin]
%
% Outputs: Beta:      The backscatter coefficient        [1/m/sr]
%          BetaTotal: The total scattering coefficient   [1/m]
%
%% Rayleigh Efficiency angular relationship
P = @(Theta) 0.7629.*(1+0.9324.*cosd(Theta).*cosd(Theta));
%% Calculating backscatter coeff
Beta = (2.938e-32).*(Press./Temp).*(1./(Lambda.^4.0117));   % Eq. (5.14)
%% Calculating the total scatter coefficient
BetaTotal = (Beta.*4.*pi)./P(180);                          % Eq. (5.15)
end
