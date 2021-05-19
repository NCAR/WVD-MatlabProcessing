% Written By: Robert Stillwell
% Written On: December 4, 2020
% Written For: National Center for Atmospheric Research

function [Spec] = ReadSystemScanData(Spectra,Const)
%
%
%
%
%
%% Modeling the etalon from known measured variables
% Determining what the spectra to be built are
[~,FN,~] = RecursiveStruct2Cell(Spectra);
% Looping over all spectra to be built
for m=1:1:length(FN)
    Spec.(FN{m,1}).Etalon = ...
        EtalonModel(Spectra.(FN{m,1}).Absorption.Lambda,Const);
end
end

function [E] = EtalonModel(Wavelength,Const)
%
% Inputs: Wavelength: An array or wavelength values for which etalon
%                     parameters should be calculated
%         Const:      A structure containing necessary universal constants
%                     
% Outputs: E:         A structure containing all of the etalon parameters
%                     needed for the WV DIAL simulation 
%
%% Defining Etalon equations
Eq.Finesse      = @(R) pi.*sqrt(R)./(1-R);
Eq.FSR          = @(L) Const.C./2./Const.IoRAir./L;
Eq.Theta        = @(L,Lambda) 4.*pi.*Const.IoRAir.*L./Lambda;
Eq.Transmission = @(R,Theta) ((1-R).^2)./(1+R.^2-2.*R.*cos(Theta));
%% Calculating etalon parameters
E.Finesse      = Eq.Finesse(Const.Etalon.Reflectivity);
E.FSR          = Eq.FSR(Const.Etalon.Length);
E.Lambda       = Wavelength;
Theta          = Eq.Theta(Const.Etalon.Length,Wavelength./1e9) - ...
                 Eq.Theta(Const.Etalon.Length,Const.Etalon.CenterWave./1e9);
E.Transmission = Eq.Transmission(Const.Etalon.Reflectivity,Theta);
end