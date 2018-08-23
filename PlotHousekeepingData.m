


function PlotHousekeepingData(PulseInfoNew,Options,Paths,PulseInfo)
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
subplot(24,1,1:4)
plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.UPS.Temperature,'r-.',...
     PulseInfoNew.TimeStamp.Merged,PulseInfoNew.WeatherStation.Temperature,'b--',...
     PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Housekeeping.Temperature,'k')
xlim([0,24]); % ylim([10,40]);
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,'Temperature [^\circC]: \color{red}UPS\color{black}, \color{blue}Weather Station\color{black}, Thermocouple', ...
     'Fontsize',FontSize)
title([Options.System, ' Housekeeping Parameters (20',Paths.Date,')'])
grid on;
set(gca,'xticklabel',{})
box on;
MaxLimits = [0,40];
YLimits = ylim;
if YLimits(1) < MaxLimits(1); YLimits(1) = MaxLimits(1); end
if YLimits(2) > MaxLimits(2); YLimits(2) = MaxLimits(2); end
ylim(YLimits);

%% Plotting other weather station data (Pressure, Relative humidity)
subplot(24,1,5:8)
[Ax, ~, ~] = plotyy(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.WeatherStation.Pressure,...
                    PulseInfoNew.TimeStamp.Merged,PulseInfoNew.WeatherStation.RelHumidity);
               
ylabel(Ax(1),'Press. [mbar]')
ylabel(Ax(2),'R.H. [%]')
xlim(Ax(1),[0,24]); xlim(Ax(2),[0,24]); %ylim(Ax(2),[0,100]);
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,'Weather Station Data','Fontsize',FontSize)
grid on;
set(gca,'xticklabel',{})

%% Plotting deviations of laser wavelength from nominal
subplot(24,1,9:12); hold on;
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
MaxLimits = [-10,10];
YLimits = ylim;
if YLimits(1) < MaxLimits(1); YLimits(1) = MaxLimits(1); end
if YLimits(2) > MaxLimits(2); YLimits(2) = MaxLimits(2); end
ylim(YLimits);

XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,'Laser \lambda Deviations [pm]','Fontsize',FontSize)
grid on;
set(gca,'xticklabel',{})
box on;


%% Plotting Laser Current
subplot(24,1,13:16); hold on;
for m=1:1:size(PulseInfoNew.Laser.Current,1)
    plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Laser.Current{m,1}.*1e3,Hue{m})
end
hold off
xlim([0,24]); 
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,'Laser Current [mA]','Fontsize',FontSize)
grid on;
set(gca,'xticklabel',{},'yaxislocation','right')
box on;

MaxLimits = [130,180];
YLimits = ylim;
if YLimits(1) < MaxLimits(1); YLimits(1) = MaxLimits(1); end
if YLimits(2) > MaxLimits(2); YLimits(2) = MaxLimits(2); end
ylim(YLimits);

%% Plotting Laser Power
subplot(24,1,17:20); hold on;
for m=1:1:size(PulseInfoNew.Laser.Power,1)
    plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Laser.Power{m,1}./1e5,Hue{m})
end
hold off
xlim([0,24]); 
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,'Laser Power x 10^5 [arb]','Fontsize',FontSize)
grid on;
set(gca,'xticklabel',{})
box on;

%% Plotting etalon temperatures
subplot(24,1,21:24); hold on;
for m=1:1:size(PulseInfoNew.Etalon.TemperatureActual,1)
    plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Etalon.TemperatureActual{m,1} - ...
                                       PulseInfoNew.Etalon.TemperatureDesired{m,1},Hue{m});
end
hold off
xlim([0,24]); 
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*TextXLoc + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*TextYLoc + YLoc(1);
text(XLoc,YLoc,'Etalon Temperature Deviations [^\circC]','Fontsize',FontSize)
xlabel('Data Time [UTC]') 
set(gca,'yaxislocation','right');
grid on;
box on;
MaxLimits = [-5,5];
YLimits = ylim;
if YLimits(1) < MaxLimits(1); YLimits(1) = MaxLimits(1); end
if YLimits(2) > MaxLimits(2); YLimits(2) = MaxLimits(2); end
ylim(YLimits);


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