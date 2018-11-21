

close all; clear; clc;

%% Loading trained PCA set
tic
cd PCAResults/
H2O    = load('H2OOnline129GHzPCA_HighResWavelength_HigherResTemp.mat');
O2On   = load('O2Online20GHzPCA.mat');
O2Off  = load('O2Offline20GHzPCA.mat');
cd ..
toc

%% Defining the temperature and pressure of interest
PressRebuild = 0.5.*ones(560,1);
TempRebuild  = linspace(230,273.15,size(PressRebuild,1))';

% % % %% Using PCA training set results to rebuild spectrum
% % % tic
% % % A = Rebuild(H2O.MeanSpectrum,H2O.PolyFitParams,H2O.PrincipleComponents, ...
% % %             PressRebuild    ,H2O.MeanP        ,H2O.SigmaP             , ...
% % %             TempRebuild     ,H2O.MeanT        ,H2O.SigmaT);
% % % B = Rebuild(O2On.MeanSpectrum,O2On.PolyFitParams,O2On.PrincipleComponents, ...
% % %             PressRebuild     ,O2On.MeanP        ,O2On.SigmaP             , ...
% % %             TempRebuild      ,O2On.MeanT        ,O2On.SigmaT);
% % % C = Rebuild(O2Off.MeanSpectrum,O2Off.PolyFitParams,O2Off.PrincipleComponents, ...
% % %             PressRebuild      ,O2Off.MeanP        ,O2Off.SigmaP             , ...
% % %             TempRebuild       ,O2Off.MeanT        ,O2Off.SigmaT);
% % % toc
% % % 
% % % %% Determining H2O absorption at exact pt(s)
% % % WavelengthsDesired = [828.188,828.3];
% % % 
% % % tic
% % % M = zeros(size(A,2),size(WavelengthsDesired,2));
% % % for m=1:1:size(PressRebuild,1)
% % %     M(m,:) = interp1(H2O.Lambda,A(:,m),WavelengthsDesired);
% % % end
% % % toc

%% Retrieving the absorption coefficients from the PCA training set
WavelengthsDesired = [828.188,828.3];
cd ..
tic
[M,A] = PCALineFitting(H2O  ,PressRebuild,TempRebuild,WavelengthsDesired);
toc
tic
[~,B] = PCALineFitting(O2On ,PressRebuild,TempRebuild,WavelengthsDesired);
toc
tic
[~,C] = PCALineFitting(O2Off,PressRebuild,TempRebuild,WavelengthsDesired);
toc
cd PCA

%% Plotting results
figure(1)
pcolor(H2O.Lambda,TempRebuild,real(log10(A'))); shading flat; colorbar
xlabel('Wavelength [nm]'); ylabel('Temperature [K]');
title('Water Vapor Absorption Coefficeint [m^2]');

figure(2)
pcolor(O2On.Lambda,TempRebuild,real(log10(B'))); shading flat; colorbar
xlabel('Wavelength [nm]'); ylabel('Temperature [K]');
title('O_2 Online Absorption Coefficeint [m^2]')

figure(3)
pcolor(O2Off.Lambda,TempRebuild,real(log10(C'))); shading flat; colorbar
xlabel('Wavelength [nm]'); ylabel('Temperature [K]');
title('O_2 Offline Absorption Coefficeint [m^2]')
