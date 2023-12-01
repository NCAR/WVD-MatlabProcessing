% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created October, 2023

function [SigmaR] = KlettInversion(Altitudes,BSCounts,Cal,GuessK,GuessSigma,GuessLR,Lambda,Pressure,Temperature)
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
%% Defining function handles
Atm2Millibar = @(X) X.*1013.25;
%% Calculating the overlap correction for MPD
Overlap = interp1(Cal.Overlap.Range,Cal.Overlap.Value,Altitudes,'linear','extrap');
%% Calculating range corrected logrithmic background subtracted counts
Z = real(log(Altitudes.*Altitudes.*BSCounts./Overlap));
% % Checking to make sure the first bin is not seeded with a zero
% BackScan = 5;   % Checking for the last ____ nan values
% for m=1:1:size(Z,2)
%     % Checking for highest valid measurements
%     Temp = find(isnan(Z(:,m)) == 0,BackScan,'last');
%     % Checking if there are enough measurements 
%     if size(Temp,1) ~= BackScan
%         Index(m,:) = [Temp;nan.*zeros(BackScan - size(Temp,1),1)];
%     else
%         Index(m,:) = Temp;
%     end
%     if size(Temp,1) ~= 0
%         % Estimating the top of the profile
%         for n=BackScan-1:-1:1
%             GuessTop(m) = Index(m,n);
%             if Index(m,n+1) - 1 ==  Index(m,n)
%                 break
%             end
%         end
%         % Filling the upper values with an estimated background
%         Z(isnan(Z(:,m)),m) = min(Z(1:GuessTop(m),m));
%     end
% end

%% Numerically integrating for the denominator term of the Klett Inversion
Integral = ExtinctionIntegral(Altitudes,Z,mean(Z(end-10:end,:)),GuessK);
%% Pre-allocating memory 
SigmaR = GuessK.*ones(size(BSCounts));
%% Iterating the Klett algorithm down for all time stamps
B = mean(Z(end-10:end,:)); 
for m=size(Altitudes,1)-1:-1:1
    SigmaR(m,:) = (exp((Z(m,:)-Z(end,:))./GuessK))./ ...
                  ((1./GuessSigma)+(2./GuessK).*Integral(m,:));
end
%% Removing points with very low SNR
SigmaR(isnan(BSCounts)) = nan;
%% Calculating the molecular contribution 
SigmaM = MolecularExtinction(Lambda,Atm2Millibar(Pressure),Temperature);
%% Calcualting the Backscatter ratio
SigmaR = SigmaR./SigmaM./GuessLR;
end

% This function numerically integrates the total scattering efficiency over
% the path to understand the total extinction caused by scattering
function [Integral] = ExtinctionIntegral (Alt,S,Sm,K)
%
% Inputs: Alt:  An array of altitudes in meters above sea level
%         S:    Range corrected counts profiles 
%         Sm:   Range corrected counts at the top of the profile
%         K:    Guess for the value of K in the constituative relationship
%               Beta = const * sigma ^ K
%         
% Outputs: Transmission: The transmission due to all scattering as a 
%                        function of height
% 
%% Resizing arrays of interest
RangeRes = Alt(2,1) - Alt(1,1);
% Sm = repmat(Sm,size(Alt,1),1);
%% Calculating extinction integral
Integral = zeros(size(Alt,1),size(S,2));
for m=size(Alt,1)-1:-1:1
    Integral(m,:) = trapz(exp(S(m:end,:)-Sm)./K,1).*RangeRes;
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
% BetaM    = Beta.*1e3;  % Converting to [km^-1]
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
