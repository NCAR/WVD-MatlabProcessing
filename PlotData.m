% Written by: Robert Stillwell
% Written for: National Center For Atmospheric Research
% This function calls the subfunctions that plot the MPD housekeeping and
% data figures. 
% Modification info: Created: November 13, 2018

function PlotData(Altitude,Counts,DataProducts,Map,Options,Paths,PulseInfo,PulseInfoNew,SurfaceWeather)
%
% Inputs: Altitude:       Structure containing all altitude arrays 
%         Counts:         Structure containing the observed photon counts
%                         and all arrays at verious stages of data 
%                         processing 
%         DataProducts:   A
%         Map:            
%         Options:        Structure containing all of the user defined 
%                         processing options
%         Paths:          A named structure containing all of the file path
%                         information needed for finding raw data and 
%                         saving processed data and figures
%         PulseInfo:      A
%         PulseInfoNew:   A
%         SurfaceWeather: Structure containing all time series information
%                         from the MPD weather stations
%
% Outputs: none
%
%% Loading colormap           
Plotting.ColorMap = importdata([Paths.Colormap,'/NCAR_C_Map.mat']);

%% 
fprintf('Plotting Data\n')
Plotting.ScreenSize = get(0,'ScreenSize');
Plotting.FontSize   = 14;
Plotting.PlotSize1  = [Plotting.ScreenSize(4)/1.5 Plotting.ScreenSize(4)/10 Plotting.ScreenSize(3)/1.5 Plotting.ScreenSize(4)/3];
Plotting.PlotSize2  = [1 Plotting.ScreenSize(4)/2 Plotting.ScreenSize(3)/2 Plotting.ScreenSize(4)/2];
Plotting.x          = (PulseInfo.DataTimeDateNumFormat)';
Plotting.xdata      = linspace(fix(min(PulseInfo.DataTimeDateNumFormat)),...
                              ceil(max(PulseInfo.DataTimeDateNumFormat)), 25);
Plotting.y          = (Altitude.RangeOriginal./1e3);

%%
PlottingMainPlots(Counts,PulseInfo,DataProducts,Options,Paths,Plotting,SurfaceWeather,Map)
 
cd(Paths.Code) % point back to original directory

PlotHousekeepingData(DataProducts,PulseInfoNew,Options,Paths,Plotting,PulseInfo)

end



function PlottingMainPlots(Counts,PulseInfo,DataProducts,Options,Paths,Plotting,SurfaceWeather,Map)
%
%
%
%
%
%
%
%% Final information needed before plotting
date = datestr(nanmean(PulseInfo.DataTimeDateNumFormat), 'dd mmm yyyy');
%% Plotting relative backscatter
figure1 = figure('visible', 'off','Position',[Plotting.ScreenSize(4)/1.5 Plotting.ScreenSize(4)/10 Plotting.ScreenSize(3)/1.5 Plotting.ScreenSize(4)/1.5]);
set(figure1, 'visible', 'off', 'PaperUnits', 'points', 'PaperPosition', [0 0 1280 800]);
if Options.flag.plot_data == 1
    set(figure1, 'visible', 'on');
end
subplot1=subplot(2,1,1,'Parent',figure1);
box(subplot1,'on');
set(gcf,'renderer','zbuffer');
Z = double(log10((real(Counts.RelativeBackscatter{Map.Offline,1}'))));
h = pcolor(Plotting.x,Plotting.y,Z);
set(h, 'EdgeColor', 'none');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.005; 0.0025]);
set(gca, 'XTick',  Plotting.xdata)
colorbar('EastOutside');
axis([fix(min(PulseInfo.DataTimeDateNumFormat)) fix(min(PulseInfo.DataTimeDateNumFormat))+1 0 12])
caxis([1 6]);
datetick('x','HH','keeplimits', 'keepticks');
colormap(Plotting.ColorMap)
%shading interp
hh = title({[date,'  Relative Backscatter (C/ns km^2)']},'fontweight','b','fontsize',Plotting.FontSize);
P_t = get(hh, 'Position');
set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
xlabel('Time (UTC)','fontweight','b','fontsize',Plotting.FontSize);
ylabel('Height (km, AGL)','fontweight','b','fontsize',Plotting.FontSize);
set(gca,'Fontsize',Plotting.FontSize,'Fontweight','b');

%% Plot water vapor in g/m^3
subplot1=subplot(2,1,2,'Parent',figure1);
box(subplot1,'on'); %(number density in mol/cm3)(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
set(gcf,'renderer','zbuffer');
Z = double(real(DataProducts.N_avg'.*1e6./6.022E23.*18.015));
h = pcolor(Plotting.x,Plotting.y,Z);
set(h, 'EdgeColor', 'none');
set(gca, 'XTick',  Plotting.xdata)
set(gca,'TickDir','out');
set(gca,'TickLength',[0.005; 0.0025]);
colorbar('EastOutside');
axis([fix(min(PulseInfo.DataTimeDateNumFormat)) fix(min(PulseInfo.DataTimeDateNumFormat))+1 0 6])
if strcmp(Options.System,'DIAL01')
caxis([0,20]);
else
caxis([0 10]);
end
datetick('x','HH','keeplimits', 'keepticks');
colormap(Plotting.ColorMap)
%shading interp
hh = title({[date,'  Water Vapor (g/m^{3})']},'fontweight','b','fontsize',Plotting.FontSize);
P_t = get(hh, 'Position');
set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
xlabel('Time (UTC)','fontweight','b','fontsize',Plotting.FontSize);
ylabel('Height (km, AGL)','fontweight','b','fontsize',Plotting.FontSize);
set(gca,'Fontsize',Plotting.FontSize,'Fontweight','b');

%% Saving the quickloook plot
if Options.flag.save_quicklook == 1
    fprintf('Making Quicklook and Uploading to Field Catalog\n')
    cd(Paths.Figures) % point to the directory where data is stored
    date=datestr(nanmean(PulseInfo.DataTimeDateNumFormat), 'yyyymmdd');
    % save the image as a PNG to the local data folder
    name=strcat('lidar.NCAR-WV-',Options.System,'_',Options.Location,'.20',Paths.Date, '0000.Backscatter_WV.png');
    print(figure1, name, '-dpng', '-r300') % set the resolution as 300 dpi
    if Options.flag.save_catalog == 1 % upload figure to the field catalog
        test=ftp('catalog.eol.ucar.edu', 'anonymous', 'spuler@ucar.edu');
        cd(test, Paths.Catalog);
        mput(test, name);
        cd(test);
        dir(test,'lidar*')
        close(test);
    end
    cd(Paths.Code)
end

    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data2Save.WavelengthX   = PulseInfo.DataTimeDateNumFormat;
% Data2Save.WavelengthOn  = PulseInfo.Lambda{1,1};
% Data2Save.WavelengthOff = PulseInfo.Lambda{2,1};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CurrentPath = pwd;
% MeasurementsTotal = size(Counts.Raw{1,1},1);
% cd(Paths.SaveData)
% ExtraSave = [Paths.FigureType,'NewDataProcessing20',num2str(Paths.FolderDate),'.mat'];
% save(ExtraSave,'Data2Save','MeasurementsTotal');
% cd(CurrentPath)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
end

function PlotHousekeepingData(DataProducts,PulseInfoNew,Options,Paths,Plotting,PulseInfo)
%
%
%
%
%
%
%% Plotting constatns
Hue = {'k','r','b','g','c','m'};
FontSize = 14;
TextXLoc = 0.02;
TextYLoc = 0.85;
screen_size = get(0,'ScreenSize');
sh = screen_size(4);                 % Screen Height
fw = 1000;                           % Figure Width
fh = 800;                            % Figure Height
mh = 90;                             % Estimated Menu Height

%% Making figure
FigureHandle = figure;
set(gcf,'Color',[1 1 1],'Position',[10 (sh-fh-mh)  fw fh]);

%% Plotting Temperatures (thermocouples, UPS, Weather station)
subplot(28,1,1:4)
plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.UPS.Temperature,'r-.',...
     PulseInfoNew.TimeStamp.Merged,PulseInfoNew.WeatherStation.Temperature,'b--',...
     PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Housekeeping.Temperature,'k')
xlim([0,24]); 
title([Options.System, ' Housekeeping Parameters (20',Paths.Date,')'])
grid on; box on;
set(gca,'xticklabel',{})
MaxLimits = [0,40];
YLimits = ylim;
if YLimits(1) < MaxLimits(1); YLimits(1) = MaxLimits(1); end
if YLimits(2) > MaxLimits(2); YLimits(2) = MaxLimits(2); end
ylim(YLimits);
AddPlotText(TextXLoc,TextYLoc,'Temperature [^\circC]: \color{red}UPS\color{black}, \color{blue}Weather Station\color{black}, Thermocouple',FontSize)


%% Plotting other weather station data (Pressure, Relative humidity)
subplot(28,1,5:8)
[Ax, ~, ~] = plotyy(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.WeatherStation.Pressure,...
                    PulseInfoNew.TimeStamp.Merged,PulseInfoNew.WeatherStation.RelativeHumidity);
ylabel(Ax(1),'Press. [mbar]'); ylabel(Ax(2),'R.H. [%]')
xlim(Ax(1),[0,24]); xlim(Ax(2),[0,24]); 
grid on; box on;
set(gca,'xticklabel',{})
AddPlotText(TextXLoc,TextYLoc,'Weather Station Data',FontSize)

%% Plotting deviations of laser wavelength from nominal
subplot(28,1,9:12); hold on;
WorstSigma = 0;
for m=1:1:size(PulseInfoNew.Laser.WavelengthActual,1)
    plot(PulseInfoNew.TimeStamp.Merged,1000.*(PulseInfoNew.Laser.WavelengthActual{m,1} - ...
                                              PulseInfoNew.Laser.WavelengthDesired{m,1}),Hue{m});
    BoundsTemp = nanstd(1000.*(PulseInfoNew.Laser.WavelengthActual{m,1} - ...
                                              PulseInfoNew.Laser.WavelengthDesired{m,1}));
    if abs(BoundsTemp) > abs(WorstSigma)
        WorstSigma = BoundsTemp;
    end
end
hold off
xlim([0,24]); 
MatlabYLimits = ylim; MatlabYLimits = max(abs(MatlabYLimits));
if MatlabYLimits > abs(6.*WorstSigma)
    ylim([-abs(6.*WorstSigma),abs(6.*WorstSigma)])
else
    ylim([-MatlabYLimits,MatlabYLimits])
end
% Checking if the matlab bounds o6 6 sigma are better
MaxLimits = [-10,10]; MinLimits = [-0.2,0.2];
YLimits = ylim;
if YLimits(1) < MaxLimits(1); YLimits(1) = MaxLimits(1); end
if YLimits(2) > MaxLimits(2); YLimits(2) = MaxLimits(2); end
ylim(YLimits);
if YLimits(1) > MinLimits(1); YLimits(1) = MinLimits(1); end
if YLimits(2) < MaxLimits(2); YLimits(2) = MinLimits(2); end
ylim(YLimits);
grid on; box on;
set(gca,'xticklabel',{})
AddPlotText(TextXLoc,TextYLoc,'Laser \lambda Deviations [pm]',FontSize)


%% Plotting Laser Current
subplot(28,1,13:16); hold on;
for m=1:1:size(PulseInfoNew.Laser.Current,1)
    plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Laser.Current{m,1}.*1e3,Hue{m})
end
hold off
xlim([0,24]); 
grid on;
set(gca,'xticklabel',{},'yaxislocation','right')
box on;

MaxLimits = [110,180];
YLimits = ylim;
if YLimits(1) < MaxLimits(1); YLimits(1) = MaxLimits(1); end
if YLimits(2) > MaxLimits(2); YLimits(2) = MaxLimits(2); end
ylim(YLimits);
AddPlotText(TextXLoc,TextYLoc,'Laser Current [mA]',FontSize)

%% Plotting Laser Power
subplot(28,1,17:20); hold on;
for m=1:1:size(PulseInfoNew.Laser.Power,1)
    plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Laser.Power{m,1}./1e5,Hue{m})
end
hold off
xlim([0,24]); 
grid on; box on;
set(gca,'xticklabel',{})
AddPlotText(TextXLoc,TextYLoc,'Laser Power x 10^5 [arb]',FontSize)

%% Plotting etalon temperatures
subplot(28,1,21:24); hold on;
for m=1:1:size(PulseInfoNew.Etalon.TemperatureActual,1)
    plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Etalon.TemperatureActual{m,1} - ...
                                       PulseInfoNew.Etalon.TemperatureDesired{m,1},Hue{m});
end
hold off
xlim([0,24]); 
set(gca,'xticklabel',{},'yaxislocation','right');
grid on; box on;
MaxLimits = [-5,5];MinLimits = [-0.01,0.01];
YLimits = ylim;
if YLimits(1) < MaxLimits(1); YLimits(1) = MaxLimits(1); end
if YLimits(2) > MaxLimits(2); YLimits(2) = MaxLimits(2); end
ylim(YLimits);
AddPlotText(TextXLoc,TextYLoc,'Etalon Temperature Deviations [^\circC]',FontSize)
DesiredWidth = get(gca,'position');
DesiredWidth = DesiredWidth(3);

%% Plotting optical depth
subplot(28,1,25:28); hold on;
contour3((PulseInfo.DataTime - PulseInfo.DataTime(1)).*24,Plotting.y,DataProducts.OpticalDepth',...
         [0:0.5:2],'ShowText','on','LineColor','k'); 
pcolor ((PulseInfo.DataTime - PulseInfo.DataTime(1)).*24,Plotting.y,DataProducts.OpticalDepth'); 
hold off
shading flat; colorbar; 
caxis([0,2]); xlim([0,24]);ylim([0,5])
set(gca,'ytick',0:1:5)
grid on; box on;
xlabel('Data Time [UTC]') 
AddPlotText(TextXLoc,TextYLoc,'2-way Water Vapor Optical Depth',FontSize)
CurrentWindow = get(gca,'position');
set(gca,'position',[CurrentWindow(1:2),DesiredWidth,CurrentWindow(4)])

%% Saving plot
if Options.flag.save_quicklook == 1
    fprintf('Making Housekeeping Figure and Uploading to Field Catalog\n')
    cd(Paths.Figures) % point to the directory where data is stored
    date=datestr(nanmean(PulseInfo.DataTimeDateNumFormat), 'yyyymmdd');
    % save the image as a PNG to the local data folder
    name=strcat('lidar.NCAR-WV-',Options.System,'_',Options.Location,'.20',Paths.Date, '0000.Housekeeping.png');
    print(FigureHandle, name, '-dpng', '-r300') % set the resolution as 300 dpi
    if Options.flag.save_catalog == 1 % upload figure to the field catalog
        test=ftp('catalog.eol.ucar.edu', 'anonymous', 'spuler@ucar.edu');
        cd(test, Paths.Catalog);
        mput(test, name);
        cd(test);
        dir(test,'lidar*')
        close(test);
    end
    cd(Paths.Code)
end

end

function AddPlotText(TextXLoc,TextYLoc,Text,FontSize)
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,Text,'Fontsize',FontSize)
end
