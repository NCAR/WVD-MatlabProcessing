


function [Constants,Conversions] = DefineConstants


%% Defining universal constants
Constants.m   = 18.015E-3./6.022E23; % mass of a single water molecule
Constants.k_B = 1.3806488e-23;       % (J/K)
Constants.c   = 299792458;           % (m/s) (exact)

%% Defining unit conversions
Conversions.Pas2Atm = @(Pascals) Pascals.*9.86923e-6;

end