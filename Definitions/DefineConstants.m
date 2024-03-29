% Written By: Robert Stillwell
% Written On: September 29, 2020
% Written For: National Center for Atmospheric Research

function [Const] = DefineConstants
%
% Inputs: none
%
% Outputs: Const: A structure containing all of the universal and earth
%                 constants needed for temperature analysis
%
%% Universal constants
Const.Av      = 6.0221408e23;   % Avogadros number          [molec/mole]
Const.C       = 299792458;      % Speed of light            [m/s] 
Const.H       = 6.626e-34;      % Planck's constant         [Js]
Const.Kb      = 1.38065e-23;    % Boltzman's constant       [m^2 kg/s^2/K]
Const.R       = 8.3144598;      % Universal gas constant    [J/(mol*K)]
%% Constants on earth
Const.G0      = 9.80665;       % Gravitational acceleration   [m/s^2]
Const.IoRAir  = 1.00037;       % Index of refraction of air   [unitless]
Const.MO2     = 5.314e-26;     % Mass O2 molecule             [kg]
Const.MWV     = 2.9915e-26;    % Mass H2O molecule            [kg]
Const.MAir    = 4.792e-26;     % Mass of air                  [kg]
Const.MolMAir = 0.0289644;     % Molar mass of Earth's air    [kg/mol]
Const.No      = 2.47937e25;    % Loschmidt's number (@296K & 1 atm) [1/m^3]
Const.QO2     = .2095;         % O2 dry mixing ratio          [unitless]
%% Conversions
Const.Atm2Pa   = 101325;        % Conversion from atmospheres to Pascals
Const.C2K      = 273.15;        % Conversion from Celsius to Kelvin
Const.MBar2Atm = 1013.25;       % Conversion from millibar to atmospheres
%% Optics
Const.Etalon.CenterWave   = 769.7958;    % nanometers
Const.Etalon.Length       = 0.00094896;  % Meters
Const.Etalon.Reflectivity = 0.816072;    % Unitless
%% Spectroscopy
Const.Eo  = [1420.766;1635.0686].*100;      % Ground state energy          [m^-1]
Const.O2LineS = [4.86e-26;3.749e-27]./100;  % Line strength of O2 line     [m/Mol]
Const.Eo  = Const.Eo*Const.H*Const.C;       % Same converted to units of   [J]
Const.HitranTo = 296;                       % Hitran reference temperature [K]
%% Reshaping Line parameters
Const.Eo      = PushTo3d(Const.Eo);
Const.O2LineS = PushTo3d(Const.O2LineS);
end

function [ThreeD] = PushTo3d(OneD)
%
% Input: OneD: Array of line parameters as 1d array
%
% Output: ThreeD: Array reshaped such that the line parameters occupy the
%                 3rd array position. Will be used as (time,range,line)
%%
ThreeD = reshape(OneD,[1,1,length(OneD)]);
end