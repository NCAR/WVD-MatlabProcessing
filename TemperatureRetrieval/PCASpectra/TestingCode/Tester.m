
%% Setting up workspace
% close all;
clear; clc
%% Function Handles
PercentEr = @(A,B)          (A-B)./A.*100;
DLam2DNu  = @(Lam,DLam)   -299792458.*DLam./Lam./Lam;
%% Loading data
Atmo   = load('TesterAtmo.mat');
cd ..
Online   = load('O2Online20GHzPCA.mat');
OnlineRB = load('O2OnlineRB20GHzPCA.mat');
cd TestingCode
%% Defining desired wavelength spectrum
CenterLambda = 769.7958;
LambdaRange  = 0.1;
WavelengthsDesired = linspace(CenterLambda-LambdaRange/2, ...
                              CenterLambda+LambdaRange/2,499);

%% Defining random profile to check
AltIn  = floor(size(Atmo.Pressure,1).*rand);
TimeIn = floor(size(Atmo.Pressure,2).*rand);
% Reshaping arrays for ease of use
TempRebuild  = Atmo.Temperature(:);
PressRebuild = Atmo.Pressure(:);
%% Recreating spectrum with principle component analysis
tic
[CrossSection,~]   = PCALineFitting(Online,PressRebuild,TempRebuild,WavelengthsDesired);
[CrossSectionRB,~] = PCALineFitting(OnlineRB,PressRebuild,TempRebuild,WavelengthsDesired);
%% Reshaping spectra
CrossSection2   = reshape(CrossSection,[size(Atmo.Pressure,1),size(Atmo.Pressure,2),size(WavelengthsDesired,2)]);
CrossSectionRB2 = reshape(CrossSectionRB,[size(Atmo.Pressure,1),size(Atmo.Pressure,2),size(WavelengthsDesired,2)]);
toc
%% Building test spectrum
FreqSpec = DLam2DNu(CenterLambda.*1e-9,(WavelengthsDesired-CenterLambda).*1e-9);
tic
CrossSection3   = BuildHitranSpectrum(WavelengthsDesired,Atmo.Pressure(AltIn,TimeIn),Atmo.Temperature(AltIn,TimeIn),'O2');
CrossSectionRB3 = RayleighBrillouinSpecWavelength(FreqSpec,CenterLambda.*1e-9,Atmo.Pressure(AltIn,TimeIn).*101325,Atmo.Temperature(AltIn,TimeIn));
toc
%% Plotting data to make sure it works
Figures = double(findall(0, 'type', 'figure'));
if isempty(Figures)
    FigNum = 0;
else
    FigNum = max(Figures);
end

figure(FigNum+1);
subplot(2,1,1)
plot(WavelengthsDesired,CrossSection3,'k',...
    WavelengthsDesired,squeeze(CrossSection2(AltIn,TimeIn,:)),'g-.')
xlabel('Wavelength'); ylabel('Absorption Cross Section [cm^2]');
set(gca,'yscale','log')
subplot(2,1,2)
plot(WavelengthsDesired,PercentEr(CrossSection3',squeeze(CrossSection2(AltIn,TimeIn,:))),'k')
xlabel('Wavelength'); ylabel('Percent Error');
title(['Temperature: ',num2str(Atmo.Temperature(AltIn,TimeIn)),', Pressure: ',num2str(Atmo.Pressure(AltIn,TimeIn))])



figure(FigNum+2);
subplot(2,1,1)
plot(WavelengthsDesired,CrossSectionRB3,'k',...
     WavelengthsDesired,squeeze(CrossSectionRB2(AltIn,TimeIn,:)),'g-.')
xlabel('Wavelength'); ylabel('Rayleigh-Brilouin Spectrum');
subplot(2,1,2)
plot(WavelengthsDesired,PercentEr(CrossSectionRB3',squeeze(CrossSectionRB2(AltIn,TimeIn,:))),'k')
xlabel('Wavelength'); ylabel('Percent Error');
title(['Temperature: ',num2str(Atmo.Temperature(AltIn,TimeIn)),', Pressure: ',num2str(Atmo.Pressure(AltIn,TimeIn))])