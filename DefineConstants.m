% Written by: Scott Spuler
% Written for: National Center For Atmospheric Research
% This function defines universal constants and conversion function handles
% for use throughout the MPD matlab processing code
% Modification info: Created: November 13, 2018

function [Constants,Conversions] = DefineConstants
%
% Inputs: none
%
% Outputs: Constants:   A structure containing named universal constants
%          Conversions: A structure containing named function handles used
%                       to convert units throughout the MPD processing code
%
%% Defining universal constants
Constants.c   = 299792458;           % (m/s) (exact)
Constants.m   = 18.015E-3./6.022E23; % mass of a single water molecule (kg)
Constants.N_A = 6.0221415E23;        %mol^-1
Constants.k_B = 1.3806488e-23;       % (J/K)
Constants.R   = 8.31447215;   %J mol^-1 K^-1

%% Defining unit conversions
Conversions.Pas2Atm         = @(Pascals) Pascals.*9.86923e-6;
Conversions.CelciusToKelvin = @(Celcius) Celcius + 273.15;

end