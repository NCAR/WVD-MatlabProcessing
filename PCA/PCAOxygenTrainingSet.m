

%% Setting up runtime environment
clear all; close all; clc;

%% Function handles
NormalizePolyInputs   = @(T1,Tm,Sigma)    (T1-Tm)./Sigma;
DNu2DLamWithNu        = @(C,Nu,DNu)       -C.*DNu./Nu./Nu;

%% Options
Online = 'no';

%% Constants
C            = 299792458;        % Speed of light in vacuum [m/s]
CenterLam    = 770.10;           % Center wavelength to calculate/plot [nm]
FreqWindow   = 20e9;             % Full width of the wavelength window [Hz]
Point2Plot   = 1999;
Pressures    = linspace(0.01,1.2,50);
Temperatures = linspace(175,325,40);
Gas          = 'O2';

[Pressures,Temperatures] = meshgrid(Pressures,Temperatures);

if strcmp(Online,'yes')
    Type = 'Online';
else
    Type = 'Offline';
%     CenterLam   = CenterLam - DNu2DLamWithNu(C,C./(CenterLam./1e9),43e9).*1e9; 
end

%% Loading Hitran data  
LambWindow = DNu2DLamWithNu(C,C./(CenterLam./1e9),FreqWindow/2).*1e9;
% Making wavelength list
Lambda = linspace(CenterLam+LambWindow,CenterLam-LambWindow,Point2Plot);
% Loading Hitran data
cd HitranData/; 
formatSpec = '%*s%*s%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fid = fopen('O2_HITRAN2012_760_781.txt');
Hitran = textscan(fid, formatSpec,'headerLines', 14,'Delimiter',',');
Hitran{end} = [];
Hitran = cell2mat(Hitran);
cd ..
% Converting wavelength to cm^-1
Nu = (1e7)./Lambda;

%% Making absorption spectra
Counter      = 1;
TotalSpectra = size(Pressures,1)*size(Pressures,2);
Sigma = zeros(TotalSpectra,Point2Plot);
for m=1:1:size(Pressures,1)
    for n=1:1:size(Pressures,2)
        Sigma(Counter,:) = BuildAbsorptionSpectrum(Hitran,Nu,Pressures(m,n),Temperatures(m,n),Gas);
        fprintf('Producing Spectrum %0.0f of %0.0f\n',Counter,TotalSpectra)
        Counter = Counter + 1;
    end
end

%% Principle component analysis
MeanSpectrum = mean(Sigma);
A            = Sigma - repmat(MeanSpectrum,size(Sigma,1),1);
[PrincipleComponents,PrincipleWeights,V2]= svd(A');
PrincipleWeights  = PrincipleWeights*V2';

%% Creating weight contours and fitting polynomials
PC2Use = 50; 

% Reshaping the temperature and pressure matrices into a vector
TFit = reshape(Temperatures',size(Temperatures,1)*size(Temperatures,2),1);
PFit = reshape(Pressures',size(Pressures,1)*size(Pressures,2),1);

% Normalization parameters
MeanT  = mean(TFit);
SigmaT = 25;
MeanP  = mean(PFit);
SigmaP = 0.2;

% Normalizing Temperature and Pressure to make them semi-Gaussian instead
% of monopolar
TFit = NormalizePolyInputs(TFit,MeanT,SigmaT);
PFit = NormalizePolyInputs(PFit,MeanP,SigmaP);

for m=1:1:PC2Use
    % Reshaping the PCA weights into a matrix instead of a vector 
    WeightContour{m,1} = reshape(PrincipleWeights(m,:),size(Pressures,2),size(Pressures,1))';
    % General order XY fit polynomial
    PolyFitParams{m,1} = polyfitn([TFit,PFit],PrincipleWeights(m,:)',10);
end

%% Random rebuilding spectrum
PC2Use       = 25;
PressRebuild = 0.83;
TempRebuild  = 234.6;

for m=1:1:PC2Use
    % Rebuilding absorption spectrum weights with polynomial fits
    PolyWeights2(m,:) = polyvaln(PolyFitParams{m,1},[NormalizePolyInputs(TempRebuild,MeanT,SigmaT),NormalizePolyInputs(PressRebuild,MeanP,SigmaP)]);
end
% Rebuilding absorption spectrum in multiple ways
RebuildPCA   = PrincipleComponents(:,1:PC2Use)*PrincipleWeights(1:PC2Use,:);
RebuildPoly2 = PrincipleComponents(:,1:PC2Use)*PolyWeights2(1:PC2Use,:);
 

% %% Plotting
% figure(100); hold on;
% for m=1:1:PC2Use
%     plot3(Lambda,m.*ones(size(PrincipleComponents(:,m),1)),PrincipleComponents(:,m));
% end
% view(3)
% xlabel('Wavelength [nm]'); ylabel('Principle Component Number');
% 
figure(10);
plot(Lambda,BuildAbsorptionSpectrum(Hitran,Nu,PressRebuild,TempRebuild,Gas)./1e4,'k', ...
     Lambda,(MeanSpectrum + RebuildPoly2(:,1)')./1e4,'g-.')
axis tight;
xlabel('Wavelength [nm]'); ylabel('Cross Section [m^2]');
legend('Hitran','Polyfitn');
set(gca,'yscale','log')

PercentError = @(Ex,Th) abs((Ex-Th)./Th).*100; 
figure(11); 
plot(Lambda, PercentError((MeanSpectrum + RebuildPoly2(:,1)')./1e4, ...
                          BuildAbsorptionSpectrum(Hitran,Nu,PressRebuild,TempRebuild,Gas)./1e4));
xlabel('Wavelength [nm]'); ylabel('Percent Error [%]'); 
axis tight
 
%% Formatting plots
Current = pwd;cd('/Users/Bobby/Documents/MATLAB');MovePlots;cd(Current)

%% Cleaning workspace
% % clear A C Counter Current DNu2DLamWithNu Gas Hitran Lam2WN LambWindow m n 
% % clear NormalizePolyInputs PC2Use PercentErrorPFit Point2Plot PressRebuild
% % clear PFit PolyWeights2 PrincipleWeights Pressures Rebuild* Temperatures 
% % clear TempRebuild TFit TotalSpectra V2 WeightContour Sigma PercentError Nu
% % clear CenterLam

%% Saving PCA results
Name = ['O2',Type,num2str(FreqWindow./1e9),'GHzPCA'];
cd PCAResults/
save(Name,'Lambda','MeanP','MeanSpectrum','MeanT','PolyFitParams','PrincipleComponents','SigmaP','SigmaT')
cd ..