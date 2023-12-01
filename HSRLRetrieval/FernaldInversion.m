% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created October, 2023

function [SigmaR] = FernaldInversion(Altitudes,BSCounts,Cal,GuessLR,Lambda,Pressure,Temperature)
%
% Inputs: Altitudes:  Altitudes of data analysis                        [m]
%                     This is an [m x 1] array with for all measurements
%         BSCounts:   Background subtracted photon counts
%                     This is an [m x n] array with all measurements

%         GuessK:     Guess for the value of K in the constituative
%                     relationship Beta = const * sigma ^ K
%         GuessSigma: Guess for the value of sigma at the top of the layer
%                     of interest (sigma = scattering cross section)
%         GuessLR:    Guess for the lidar ratio                        [sr]
%
% Outputs: SigmaR:    Retrived values of total attenuation coefficient
%                     devided by the molecular attenuation coefficient
%                     including molecules and aerosols.          [unitless]
%
%% Defining function handles and Constants
Atm2Millibar = @(X) X.*1013.25;
Mol_LR = (8.*pi./3); % Lidar ratio for molecules
%% Calculating the overlap correction for MPD
Overlap = interp1(Cal.Overlap.Range,Cal.Overlap.Value,Altitudes,'linear','extrap');
%% Calculating range corrected logrithmic background subtracted counts
Z = real(Altitudes.*Altitudes.*BSCounts./Overlap);
%% Calculating the molecular contribution 
SigmaM = MolecularExtinction(Lambda,Atm2Millibar(Pressure),Temperature);
%% Numerically integrating terms for Fernald Inversion
Integral_A = DownwardIntegral(Altitudes,SigmaM);
Integral_B = DownwardIntegral(Altitudes,Z.*exp(-2.*(GuessLR-Mol_LR).*(-1.*SigmaM)));
%% Calculating aerosol backscatter coefficient 
SigmaR = (Z.*exp(2.*(GuessLR-Mol_LR).*Integral_A))./ ...
         ((Z(end,:)./(0+SigmaM(end,:))) + 2.*GuessLR.*Integral_B);
%% Removing points with very low SNR
SigmaR(isnan(BSCounts)) = nan;
end

% This function numerically integrates the total scattering efficiency over
% the path to understand the total extinction caused by scattering
function [Integral] = DownwardIntegral (Alt,S)
%
% Inputs: Alt:  An array of altitudes in meters above sea level
%         S:    Thing to be integrated
%         
% Outputs: Integral: Pre-integrated element
% 
%% Resizing arrays of interest
RangeRes = Alt(2,1) - Alt(1,1);
%% Calculating extinction integral
Integral = zeros(size(Alt,1),size(S,2));
for m=size(Alt,1)-1:-1:1
    Integral(m,:) = trapz(S(m:end,:),1).*RangeRes;
end
end

% Calculates molecular extinction at 532 nm for Summit in [km^-1]
function [BetaM] = MolecularExtinction(Lambda,Press,Temp)
%
% Inputs: Lambda:  Laser wavelength                   [meters]
%         Press:   Atmospheric pressure               [millibar]
%         Temp:    Atmospheric temperature            [Kelvin]
%
% Outputs: BetaM:  An array of molecular extinction values based on the
%                  modeled temperature and pressure
%
%% Calculating rayleigh backscattering cross section
[BetaM,~] = RayleighBackscatterCoeff(Lambda,Press,Temp);
end

% This function calculated the Rayleigh scattering efficiency as a function
% of the wavelength of probing radiation and the atmospheric parameters
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
