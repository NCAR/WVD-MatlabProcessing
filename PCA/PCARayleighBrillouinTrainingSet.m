
%% Setting up runtime environment
clear; close all; clc;

%% Constants
Kb                  = 1.3806504e-23;
Mair                = (1.66053886e-27)*28.013;
Viscosity           = 17.63e-6;
BulkViscosity       = Viscosity*0.73;
ThermalConductivity = 25.2e-3;

% % %% Lidar constants
% % Lambda         = 828e-9;

%% Unknown Tenti Model constants
c_int     = 1.0;
c_tr      = 3/2;
gamma_int = c_int/(c_tr+c_int);
rlx_int   = 1.5*BulkViscosity/(Viscosity*gamma_int);
eukenf    = Mair*ThermalConductivity/(Viscosity*Kb*(c_tr+c_int));

%% Atmospheric Pressure and temperature
Pressure    = 1.*1.01325e5;
Temperature = 273;

%% Calculating Tenit Model
X = linspace(-7.5,7.5,500);
Y = linspace(0.01,2.1,100);
% Pre-allocating data
sptsig6      = zeros(size(Y,2),size(X,2));
TentiRunTime = zeros(size(Y,2),1);
for m=1:1:size(Y,2)
    fprintf('Tenti model run %0.0f \n',m)
    tic
    [~,sptsig6(m,:)]=crbs6(Y(m),rlx_int,eukenf,c_int,c_tr,X);
    TentiRunTime(m) = toc;
end

%% Singular value decomposition 
MeanSpectrum = sptsig6(floor(size(sptsig6,1))/2,:);
A            = sptsig6 - repmat(MeanSpectrum,size(sptsig6,1),1);
[PrincipleComponents,PrincipleWeights,V2]    = svd(A');

PaperWeights  = PrincipleWeights*V2';

%% Polynomial fitting the weights
for m=1:1:size(PrincipleComponents,1)
    PolyCoefficients(m,:) = polyfit(Y,PaperWeights(m,:),6);
    PolyWeights(m,:)      = polyval(PolyCoefficients(m,:),Y);
end
for m=1:1:size(PrincipleComponents,1)
    PolyWeights(m,:) = polyval(PolyCoefficients(m,:),Y);
end

%% Reconstructing 
PC2Use = 6; 
ReconstructPoly = PrincipleComponents(:,1:PC2Use)*PolyWeights(1:PC2Use,:);

%% Plotting
figure(1); hold on;
for m=1:1:5
    plot3(X,m.*ones(size(PrincipleComponents(:,m),1)),PrincipleComponents(:,m));
end
xlim([-5,5]); view(3)
xlabel('X'); ylabel('Principle Component Number');

figure(2);
plot(X,sptsig6(1,:)  ,'k',X,MeanSpectrum'+ReconstructPoly(:,1),'r-.', ...
     X,sptsig6(end,:),'b',X,MeanSpectrum'+ReconstructPoly(:,end),'c-.');
title('Reconstructed Spectrum from Polynomial Fit Weights'); 
xlabel('X'); ylabel('Intensity [a.u.]')

%% Formatting plots
Current=pwd;cd('/Users/Bobby/Documents/MATLAB');MovePlots;cd(Current)

%% Saving data
% cd PCAResults/
% save('RBSpectrumPCA','MeanSpectrum','PolyCoefficients','PrincipleComponents','X','Y');
% cd ..

