% Written by: Robert Stillwell
% Written for: NCAR ASP
% Revision Info: Created 09/08/2017

function [Sigma] = BuildAbsorptionSpectrum(Hitran,Nu,Pressure,Temperature,Molecule)
%
% In
%
%
%
%
%
%
%% Loading universal constants
Constants = DefineUniversalConstants(Molecule);
%% Hitran/Voigt profile reference constants
ReferencePressure    = 1;       % Hitran reference pressure [atm]
ReferenceTemperature = 296;     % Hitran reference temperature [K]
%% Parsing up hitran information
Delta            = Hitran(:,8); % Pressure shift from HiTRAN [cm^-1 atm^-1]
LineCenterRef    = Hitran(:,1); % The frequency of the absorption line  
                                % center from Hitran [cm^-1]
LineEnergy       = Hitran(:,6); % Ground state transition energy from 
                                % Hitran 
LineStrengthRef  = Hitran(:,2); % Reference line strength from Hitran 
                                % [cm^-1/(mol*cm^-2)] [cm^-1]
LineWidthAirRef  = Hitran(:,4); % Air-broadened halfwidth at reference 
                                % temperature and pressure from Hitran 
                                % [cm^-1/atm]
N                = Hitran(:,7); % Linewidth temperature dependence factor 
                                % from HITRAN
%% Shifting the line center reference location based on pressure 
LineCenterRef    = LineCenterRef+Delta.*(Pressure./ReferencePressure);
%% Calculating parameters for the Voigt profile
AlphaD = HalfWidthDoppler(Constants,LineCenterRef,Temperature);
AlphaL = HalfWidthLorentz(LineWidthAirRef,N,Pressure,ReferencePressure, ...
                                        ReferenceTemperature, Temperature);
LStrength = LineStrength(Constants,LineEnergy,LineStrengthRef,          ...
                                        ReferenceTemperature, Temperature);

%% Calculating Voigt profiles                                  
[~,Sigma]=VoigtProfile(AlphaD,AlphaL,Nu,LineCenterRef,LStrength); 

end

function [Const] = DefineUniversalConstants(Molecule)
%
%
%
%
%%
Const.C   = 299792458;               % Speed of light [m/s] (exact)
Const.h   = 6.62607004e-34;          % Planck's constant [m^2kg/s]
if strcmp(Molecule,'O2')
    Const.M   = 31.9988E-3./6.022E23;
else
    Const.M   = 18.015E-3./6.022E23; % Mass of a single water molecule
end
Const.Kb  = 1.3806488e-23;           % Boltzmann's constant [J/K]
end

function [AlphaD] = HalfWidthDoppler(Const, Nuij, Temperature)
%
% Inputs: Const:       A structure containing needed universal constants in
%                      standard MKS units
%         Nuij:        The wavenumber of the spectral line transition of
%                      interest in vacuum with units of [cm^-1]
%         Temperature: The temperature of the transition of interest with
%                      units of [K] to match the units in the constants
%                      structure
%                      
% Outputs: AlphaD:     Doppler broadening half width of a particular 
%                      transition with units of [cm^-1]
%
%% Doppler broadening half width at half maximum
AlphaD = (Nuij./Const.C).*sqrt(2.*Const.Kb.*Temperature.*log(2)./Const.M);
end

function [AlphaL] = HalfWidthLorentz(LineWidthR, N, Pressure,RefP, RefT, Temperature)
%
% Inputs: Linewidth:   The Lorentz broadening half width of the transition
%                      at standard temperature and pressure
%         N:           Linewidth temperature dependence exponent from
%                      Hitran 
%         Pressure:    The pressure of the desired Lorentz broadeding
%                      calculation with units of [atm]
%         RefP:        The reference pressure of the tabulated Lorentz
%                      broadening halfwidth with units of [atm]
%         RefT:        The reference pressure of the tabulated Lorentz
%                      broadening halfwidth with units of [K]
%         Temperature: The temperature of the desired Lorentz broadeding
%                      calculation with units of [K]
%                      
% Outputs: AlphaL:     Lorentz broadening half width of a particular 
%                      transition with units of [cm^-1]
%
%% Lorentz broadening half width at half maximum
AlphaL = LineWidthR.*(Pressure./RefP).*(RefT./Temperature).^N;
end

function [LS] = LineStrength(Const, Energy, LS0, RefT, Temperature)
%
% Inputs: Const:       A structure containing needed universal constants in
%                      standard MKS units
%         Energy:      The energy of the selected transition with units of
%                      [cm^-1]
%         LS0:         The tabulated line strength of the transition from
%                      Hitran at standard temperature and pressure
%         RefT:        The reference pressure of the tabulated Lorentz
%                      broadening halfwidth with units of [K]
%         Temperature: The temperature of the desired Lorentz broadeding
%                      calculation with units of [K]
%
% Outputs: LS:         The line strength of the transition at the desired
%                      temperature
%
%% Calculating the line strength
LS = LS0.*((RefT./Temperature).^1.5).*exp((100.*Const.h.*Const.C./Const.Kb) ...
                                     .*Energy.*(1./RefT - 1./Temperature));
end

function [Voigt,Sigma] = VoigtProfile(AlphaD, AlphaL, Nu, Nu0, S)
%
% Inputs: AlphaD: Doppler broadening half width of a particular transition
%                 with units of [cm^-1]
%         AlphaL: Lorentz broadening half width of a particular transition
%                 with units of [cm^-1]
%         Nu:     The frequency of interest
%         Nu0:    The centerline frequency of a particular line
%         S:      Line strength of each line of interest
%         IntNum: The number of integration points to perform the Voigt
%                 subintegration 
%                 
% Outputs: Voigt: The value of the Voigt profile at a particular wave
%                 number with a particular temperature and pressure
%
%% Solving for the real portion of the Voigt profile
SigmaD = AlphaD./sqrt(2*log(2));
Voigt  = real(fadf(((Nu-Nu0) + 1j.*AlphaL)./SigmaD./sqrt(2)));%./(SigmaD.*sqrt(2.*pi));

Sigma = sum((S./AlphaD).*(sqrt(log(2)./pi)).*Voigt);

end

function FF = fadf(z)

%     This program file computes the complex error function, also known as
% the Faddeeva function. The algorithmic implementation utilizes the
% approximations based on the Fourier expansion [1, 2] and the Laplace
% continued fraction [3]. The code covers with high-accuracy the entire
% complex plain required for practical applications (see also optimized C++
% source code from the RooFit package in the work [4]).
%     The code remains highly accurate even at vanishing imaginary argument
% y -> 0, where y = Im[z]. The worst detected accuracy is ~10^-13.
%
% REFERENCES
% [1] S. M. Abrarov and B. M. Quine, Efficient algorithmic implementation
%     of the Voigt/complex error function based on exponential series
%     approximation, Appl. Math. Comput., 218 (2011) 1894-1902.
%     http://doi.org/10.1016/j.amc.2011.06.072
%
% [2] S. M. Abrarov and B. M. Quine, On the Fourier expansion method
%     for highly accurate computation of the Voigt/complex error function
%     in a rapid algorithm, arXiv:1205.1768v1 (2012).
%     http://arxiv.org/abs/1205.1768
%
% [3] W. Gautschi, Efficient computation of the complex error function,
%     SIAM J. Numer. Anal., 7 (1970) 187-198.
%     http://www.jstor.org/stable/2949591
%
% [4] T. M. Karbach, G. Raven and M. Schiller, Decay time integrals in
%     neutral meson mixing and their efficient evaluation,
%     arXiv:1407.0748v1 (2014).
%     http://arxiv.org/abs/1407.0748
%
%     The code is written by Sanjar M. Abrarov and Brendan M. Quine, York
% University, Canada, September 2014. Last modifications to the code were
% made on July 2016 (see the file 'readme.txt' for more information).

ind_neg = imag(z)<0; % if some imag(z) values are negative, then ...
z(ind_neg) = conj(z(ind_neg)); % ... bring them to the upper-half plane

FF = zeros(size(z)); % define array

ind_ext  = abs(z)>8; % external indices
ind_band = ~ind_ext & imag(z)<5*10^-3; % narrow band indices

FF(~ind_ext & ~ind_band) = ...
    fexp(z(~ind_ext & ~ind_band)); % internal area. This area is ...
    % ... the most difficult for accurate computation
FF(ind_ext) = contfr(z(ind_ext)); % external area
FF(ind_band) = smallim(z(ind_band)); % narrow band

    function FE = fexp(z,tauM,maxN) % Fourier expansion approximation, ...
        % ... see [1, 2] for more information

        if nargin == 1 % assign default paramenetrs tauM and maxN
            tauM = 12; % the margin value
            maxN = 23; % the upper limit summation integer
        end

        n = 1:maxN; % initiate an array
        aN = 2*sqrt(pi)/tauM*exp(-n.^2*pi^2/tauM^2); % Fourier coefficients

        z1 = exp(1i*tauM*z); % define first repeating array
        z2 = tauM^2*z.^2; % define second repeating array

        FE = sqrt(pi)/tauM*(1 - z1)./z2; % initiate array FE
        for n = 1:maxN
            FE = FE + (aN(n)*((-1)^n*z1 - 1)./(n^2*pi^2 - z2));
        end
        FE = 1i*tauM^2*z/sqrt(pi).*FE;
    end

    function CF = contfr(z) % the Laplace continued fraction approximation

        bN = 11; % initial integer
        bN = 1:bN;
        bN = bN/2;

        CF = bN(end)./z; % start computing from the last bN
        for k = 1:length(bN) - 1
            CF = bN(end-k)./(z - CF);
        end
        CF = 1i/sqrt(pi)./(z - CF);
    end

    function SIm = smallim(z) % approximation at small imag(z)

        ind_0 = abs(real(z))<5*1e-3;

        % If |Re[z]| < 5*1e-3, then:
        SIm(ind_0) = small_z(z(ind_0));

        x = real(z); % define the repeating array
        ind_poles = false(size(x)); % initiate the array of indices

        k = 1; % the counter
        while k <= 23
            % These indices are to avoid the poles that can strongly ...
            % ... deteriorate the accuracy in computation
            ind_poles = ind_poles | abs(x - k*pi/12)<1e-4;
            k = k + 1; % just to increment the counter
        end

        % Else if |Re[z]| >= 5*1e-3, then:
        SIm(~ind_0 & ~ind_poles) = narr_band(z(~ind_0 & ~ind_poles));
        % -----------------------------------------------------------------
        % Note that the margin value tauM in the line below is taken as ...
        % 12.1 instead of the default value 12. This excludes all poles ...
        % in computation even if Im[z] -> 0.
        SIm(~ind_0 & ind_poles) = narr_band(z(~ind_0 & ind_poles),12.1,23);
        % -----------------------------------------------------------------

        function SZ = small_z(z)

            % This equation improves accuracy near the origin. It is ...
            % obtained by the Maclaurin series expansion.

            % Define the repeating arrays
            zP2=z.^2;
            zP4=zP2.^2;
            zP6=zP2.*zP4;

            SZ = (((6 - 6*zP2 + 3*zP4 - zP6).*(15*sqrt(pi) + ...
                1i*z.*(30 + 10*zP2 + 3*zP4)))/(90*sqrt(pi)));
        end

        function NB = narr_band(z,tauM,maxN) % the narrow band

            % This is just an alternative representation of the ...
            % equation (14) from [1].

            if nargin == 1 % define default parameters
                    tauM = 12; % the margin value
                    maxN = 23; % the upper limit summation integer
            end

            n = 1:maxN; % initiate an array
            aN = 2*sqrt(pi)/tauM*exp(-n.^2*pi^2/tauM^2); % The Fourier ...
                                         % ... expansion coeffisients

            z1 = cos(tauM*z); % define first repeating array
            z2 = tauM^2*z.^2; % define second repeating array

            NB = 0; % initiate the array NB
            for n = 1:maxN
                NB = NB + (aN(n)*((-1)^n*z1 - 1)./(n^2*pi^2 - z2));
            end
            NB = exp(-z.^2) - 1i*((z1 - 1)./(tauM*z) - ...
                tauM^2*z/sqrt(pi).*NB);
        end
    end

% Convert for negative imag(z) values
FF(ind_neg) = conj(2*exp(-z(ind_neg).^2) - FF(ind_neg));
end