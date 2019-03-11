 

close all; clc; clear
tic

%% Function handles
Atm2Pascal   = @(ATM)   ATM.*1.01325e5;
PercentError = @(Ex,Th) abs((Ex-Th)./Th).*100; 

%% Loading the training set
cd PCAResults/
load('RayleighBrillouinXY.mat'); 
cd ..

%% Constants
RebuildTemp   = linspace(205,300,288);
RebuildPress  = linspace(Atm2Pascal(0.01),Atm2Pascal(1),280);

%% Random rebuilding spectrum
PC2Rebuild   = 10;
RebuildFreq  = linspace(-5e9,5e9,250);
CenterLam    = 770.1085e-9;     % Center wavelength to calculate/plot [nm]
 
%% Making a mesh of temperature and pressure
[RebuildPress,RebuildTemp] = meshgrid(RebuildPress,RebuildTemp);

%% Calculating rebuilt Tenti parameters
[X,Y,~,~]   = CalculateTentiParametersNDim(RebuildPress,RebuildTemp,RebuildFreq,CenterLam,Const);
X           = shiftdim(X,size(size(X),2)-1);
RebuildPoly = RebuildNDim(MeanSpectrum,PolyFitParams,PrincipleComponents,Y);
toc

%% Interpolating spectra back to desired frequency grid
tic
RebuildPolyQ = zeros(size(X)); % Pre-allocating data storage array
for m=1:1:size(X,3)       % Looping over the columns
    for n=1:1:size(X,2)   % Looping over the rows
        % Interpolating retrieved spectrum to frequency grid of interest
        RebuildPolyQ(:,n,m) = interp1(XLimits',squeeze(RebuildPoly(:,n,m)),squeeze(X(:,n,m)));      
    end
end
% Normalizing all spectra
RebuildPolyQ = RebuildPolyQ./trapz(RebuildFreq./1e9,RebuildPolyQ);

clear m n
toc

%% Loading potassium and etalon data
tic
% Loading data
cd Data
A = load('Temp115_TemperatureScan_20181012_211836_31.541000_0.100000_101_20.750000_1.000000_1_0000.txt');
cd ..
WavelengthScan = A(2:end,5);
SeedScan       = A(2:end,6); 
IntensityScan  = A(2:end,7);
% Calculating frequency and intensity
SeedScan  = 10.^(SeedScan./10);
DLam      = WavelengthScan - CenterLam.*1e9;
DNu       = -Const.C.*(DLam./1e9)./((CenterLam).^2);
Potassium = IntensityScan./SeedScan;
Potassium = Potassium./max(Potassium);
Potassium(1:2) = 1; Potassium(end-1:end) = 1;
Potassium = interp1(DNu,Potassium,RebuildFreq,'linear','extrap')';
clear A DLam DNu IntensityScan SeedScan WavelengthScan
 
% Loading data
Bounds = 100:200;
cd Data
A = load('TemperatureScan_20181115_002226_28.532000_5.000000_200_25.000000_0.500000_10_0002.txt');
cd ..
WavelengthScan = A(Bounds,5);
SeedScan       = A(Bounds,6); 
IntensityScan  = A(Bounds,7);
% Calculating frequency and intensity
SeedScan  = 10.^(SeedScan./10);
DLam      = WavelengthScan - CenterLam.*1e9;
DNu       = -Const.C.*(DLam./1e9)./((CenterLam).^2);
Etalon    = IntensityScan./SeedScan;
Etalon    = interp1(DNu,Etalon./max(Etalon),RebuildFreq,'spline','extrap')';
clear A DLam DNu IntensityScan SeedScan WavelengthScan
toc

%% Determining the cell transmission efficiency 
tic
CombinedEfficiency   = squeeze(trapz(RebuildFreq./1e9,RebuildPolyQ.*Etalon));
MolecularEfficiency  = squeeze(trapz(RebuildFreq./1e9,RebuildPolyQ.*Potassium.*Etalon));
ExtinctionEfficiency = interp1(RebuildFreq,Potassium,0,'linear')';
toc
 

%% Plotting
figure(1);
plot(RebuildFreq./1e9,RebuildPolyQ(:,1,1),'k',RebuildFreq./1e9,Potassium,'b', RebuildFreq./1e9,Etalon,'r')
xlabel('Frequency [GHz]'); ylabel('Intensity')
legend('RB Spectrum','K Transmission','Etalon Transmission'); 
 
figure(2);
plot(RebuildFreq./1e9,RebuildPolyQ(:,1,1).*Potassium,'k')
xlabel('Frequency [GHz]'); ylabel('Intensity')

figure(3);
pcolor(RebuildPress,RebuildTemp-273.15,MolecularEfficiency); shading flat; colorbar
xlabel('Pressure [Pa]'); ylabel('Temperature [^\circC]')
title('Molecular Efficiency')

figure(4);
pcolor(RebuildPress,RebuildTemp-273.15,CombinedEfficiency); shading flat; colorbar
xlabel('Pressure [Pa]'); ylabel('Temperature [^\circC]')
title('Combined Efficiency')

