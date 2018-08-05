


function PlotHousekeepingData(PulseInfoNew,Options,Paths)


Hue = {'k','r','b','g','c','m'};

figure(10)
%% Plotting Temperatures (thermocouples, UPS, Weather station)
subplot(24,1,1:4)
plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.UPS.Temperature,'r-.',...
     PulseInfoNew.TimeStamp.Merged,PulseInfoNew.WeatherStation.Temperature,'b--',...
     PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Housekeeping.Temperature,'k')
ylim([10,40]); xlim([0,24]);
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*0.02 + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*0.90 + YLoc(1);
text(XLoc,YLoc,'Temperature [^\circC]: \color{red}UPS\color{black}, \color{blue}Weather Station\color{black}, Thermocouple')
title([Options.System, ' Housekeeping Parameters (20',Paths.Date,')'])
grid on;
set(gca,'xticklabel',{})
box on;

%% Plotting other weather station data (Pressure, Relative humidity)
subplot(24,1,5:8)
[Ax, ~, ~] = plotyy(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.WeatherStation.Pressure,...
                    PulseInfoNew.TimeStamp.Merged,PulseInfoNew.WeatherStation.RelHumidity);
               
ylabel(Ax(1),'Press. [mbar]')
ylabel(Ax(2),'R.H. [%]')
xlim(Ax(1),[0,24]); xlim(Ax(2),[0,24]); ylim(Ax(2),[0,100]);
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*0.02 + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*0.90 + YLoc(1);
text(XLoc,YLoc,'Weather Station Data')
grid on;
set(gca,'xticklabel',{})

%% Plotting deviations of laser wavelength from nominal
subplot(24,1,9:12); hold on;
for m=1:1:size(PulseInfoNew.Laser.WavelengthActual,1)
    plot(PulseInfoNew.TimeStamp.Merged,1000.*(PulseInfoNew.Laser.WavelengthActual{m,1} - ...
                                              PulseInfoNew.Laser.WavelengthDesired{m,1}),Hue{m});
end
hold off
xlim([0,24]); 
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*0.02 + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*0.90 + YLoc(1);
text(XLoc,YLoc,'Laser \lambda Deviations [pm]')
grid on;
set(gca,'xticklabel',{})
box on;


%% Plotting Laser Current
subplot(24,1,13:16); hold on;
for m=1:1:size(PulseInfoNew.Laser.Current,1)
    plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Laser.Current{m,1},Hue{m})
end
hold off
xlim([0,24]); 
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*0.02 + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*0.90 + YLoc(1);
text(XLoc,YLoc,'Laser Current [A]')
grid on;
set(gca,'xticklabel',{},'yaxislocation','right')
box on;

%% Plotting Laser Power
subplot(24,1,17:20); hold on;
for m=1:1:size(PulseInfoNew.Laser.Power,1)
    plot(PulseInfoNew.TimeStamp.Merged,PulseInfoNew.Laser.Power{m,1},Hue{m})
end
hold off
xlim([0,24]); 
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*0.02 + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*0.90 + YLoc(1);
text(XLoc,YLoc,'Laser Power [arb]')
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
XLoc = xlim; XLoc = (XLoc(2)-XLoc(1)).*0.02 + XLoc(1);
YLoc = ylim; YLoc = (YLoc(2)-YLoc(1)).*0.90 + YLoc(1);
text(XLoc,YLoc,'Etalon Temperature Deviations [^\circC]')
xlabel('Data Time [UTC]') 
set(gca,'yaxislocation','right');
grid on;
box on;


end