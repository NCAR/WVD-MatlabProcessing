

%% Setting up runtime environment
clear all; close all; clc;

%% Function handles
NormalizePolyInputs   = @(T1,Tm,Sigma)    (T1-Tm)./Sigma;
DNu2DLamWithNu        = @(C,Nu,DNu)       -C.*DNu./Nu./Nu;
Atm2Pascal            = @(ATM)            ATM.*1.01325e5;

%% Defining the universal constants
Const.AMU           = 1.66053886e-27;       % Atomic mass unit [kg]
Const.C             = 299792458;            % Speed of light in vacuum [m/s]
Const.Kb            = 1.3806504e-23;        % Boltzmann's Constant [m^2*kg/s^2/K]
Const.MAir          = 28.013.*Const.AMU;    % Mass of air [kg]
Const.Theta         = pi;                   % Scattering angle of pi (Backscattering)
Const.Viscosity     = 17.63e-6;             % Shear viscosity of air [Pa*s]
Const.ViscosityBulk = Const.Viscosity*0.73; % 
Const.TConductivity = 25.2e-3;              % Thermal conductivity 

%% Operational Constants
CenterLam    = 770.1085e-9;     % Center wavelength to calculate/plot [nm]
FreqWindow   = 10e9;            % Full width of the wavelength window [Hz]
Point2Plot   = 300;
YLimits      = linspace(0,1.2,100);
XLimits      = linspace(-8,8,Point2Plot);
Gas          = 'Air';

%% Determining the wavelength and frequency spectra 
% Determining the wavelength window
LambWindow = DNu2DLamWithNu(Const.C,Const.C./CenterLam,FreqWindow/2);
% Making wavelength list
Lambda     = linspace(CenterLam+LambWindow,CenterLam-LambWindow,Point2Plot);

%% Making Rayleigh Brillouin spectra
Counter      = 1;
TotalSpectra = size(YLimits,1)*size(YLimits,2);
Sigma = zeros(TotalSpectra,Point2Plot);
for m=1:1:size(YLimits,2)
        Sigma(Counter,:) = BuildRBSpectrumXY(Const,XLimits,YLimits(m));
        fprintf('Producing Spectrum %0.0f of %0.0f\n',Counter,TotalSpectra)
        Counter = Counter + 1;
end

%% Principle component analysis
MeanSpectrum = mean(Sigma);
A            = Sigma - repmat(MeanSpectrum,size(Sigma,1),1);
[PrincipleComponents,PrincipleWeights,V]= svd(A');
PrincipleWeights  = PrincipleWeights*V';

%% Creating weight contours and fitting polynomials
PC2Use = 15; 

for m=1:1:PC2Use
    % Fit polynomial
    PolyFitParams{m,1} = polyfit(YLimits,PrincipleWeights(m,:),10);
end


%% Random rebuilding spectrum
PC2Rebuild   = 10;
PressRebuild = Atm2Pascal(0.2);
TempRebuild  = 230;
FreqRebuild  = linspace(-5e9,5e9,500);
CenterLam    = 770.1085e-9;     % Center wavelength to calculate/plot [nm]

% Building spectrum and interpolating to the desired X grid
[X,Y,~,~]    = CalculateTentiParametersNDim(PressRebuild,TempRebuild,FreqRebuild,CenterLam,Const);
RebuildPoly2 = RebuildNDim(MeanSpectrum,PolyFitParams,PrincipleComponents,Y);
RebuildPoly2 = interp1(XLimits,RebuildPoly2',X);
RebuildPoly2 = RebuildPoly2./trapz(FreqRebuild./1e9,RebuildPoly2);

% Running the Tenit model to check
Truth = BuildRBSpectrum(Const,FreqRebuild,Lambda,CenterLam,PressRebuild,TempRebuild);
Truth = Truth./trapz(FreqRebuild./1e9,Truth);


%% Plotting
% Checking fit contours
for m=1:1:6
   figure(1); subplot(6,1,m)
   plot(YLimits,PrincipleWeights(m,:),'k.',YLimits,polyval(PolyFitParams{m,1},YLimits),'b') 
   if m == 1
       legend('Training Set','Fit Polynomials')
   elseif m==6
       xlabel('Y Parameter'); 
   end
   ylabel('Weight')
end

figure(100); hold on;
for m=1:1:5 %PC2Use
    plot3(Lambda,m.*ones(size(PrincipleComponents(:,m),1)),PrincipleComponents(:,m));
end
view(3)
xlabel('X'); ylabel('Principle Component Number');

figure(10);
plot(FreqRebuild,Truth,'k',FreqRebuild,RebuildPoly2,'g-.')
axis tight;
xlabel('Wavelength [nm]'); ylabel('Cross Section [m^2]');
legend('Tenti 6','Polyfitn')
 
PercentError = @(Ex,Th) abs((Ex-Th)./Th).*100; 
figure(11); 
plot(FreqRebuild, PercentError(RebuildPoly2,Truth));

%% Cleaning up workspace
clear A Atm2Pascal CenterLam Counter DNu2DLamWithNu FreqRebuild FreqWindow
clear Gas Lambda LambWindow m NormalizePolyInputs PC2Rebuild PC2Use PercentError
clear Point2Plot PressRebuild PrincipleWeights RebuildPoly2 Sigma TempRebuild
clear TotalSpectra Truth  X Y 

%% Saving data 
cd PCAResults/
save('RayleighBrillouinXY','Const','MeanSpectrum','PolyFitParams','PrincipleComponents','V','XLimits','YLimits')
cd ..

