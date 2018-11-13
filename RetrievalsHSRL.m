% Written by: Robert Stillwell
% Written for: National Center For Atmospheric Research
% This function processes all the High Spectral Resolution lidar retrievals 
% from the MPD system. It also filters and averages the data.
% Modification info: Created: November 13, 2018

function [] = RetrievalsHSRL(Altitude, Capabilities, DataProducts, Map, Options, PulseInfo, SurfaceWeather)
%
%
%
%
%
%
% 
%
%% Finding universal constants
[Constants,Conversions] = DefineConstants;

%% Defining function handles (backscatter coefficient in m^-1 sr^-1)
Beta = @(Pressure,Temperature,Wavelength) 5.45.*10.^-32.*(550./Wavelength).^4.*Pressure./(Temperature.*Conversions.Pas2Atm(Constants.k_B)); %backscatter coefficient in m^-1 sr^-1

%% Molecular backscatter 
fprintf('      HSRL Retrieval: Calculating Molecular Profile\n')
% Calculate temperature and pressure profile based on surface measurement and 
% asuming a standard lapse rate (-6.5 deg/km) for the entire troposphere
if Options.flag.WS == 1
    T0 = SurfaceWeather.Temperature+273.15; % Converting from C to Kelvin
    P0 = SurfaceWeather.Pressure;
    % Adding median data to empty weather station data
    if sum(isnan(T0)) > 0
       T0(isnan(T0)) = nanmean(T0);
       P0(isnan(P0)) = nanmean(P0);
    end
else
	T0 = 290;  % Surface temperature
    P0 = 0.83; % Surface pressure in Boulder
end
% Standard moist adiabatic lapse rate (-6.5 deg/km) 
DataProducts.TemperatureProfile = T0-0.0065.*Altitude.RangeOriginal; 
% Hydrostatic equation and ideal gas law
DataProducts.PresssureProfile   = P0.*(T0./DataProducts.TemperatureProfile).^-5.5;   
% Calculating the molecular backscatter profile
DataProducts.BetaMProfile   = Beta(DataProducts.PresssureProfile,   ...
                                   DataProducts.TemperatureProfile, ...
                                   PulseInfo.LambdaNearest{Map.Molecular,1});
clear Beta P0 T0

%% Aerosol backscatter coefficient
if Capabilities.O2HSRL == 1
    fprintf('      HSRL Retrieval: Calculating Potassium HSRL Backscatter Coefficient\n')
else
    fprintf('      HSRL Retrieval: Calculating Standard HSRL Backscatter Coefficient\n')
end

% % % % %% Calculate HSRL parameters
% % % % R_size                = 150;
% % % % % Calculating the HSRL retrievals
% % % % [DataProducts.AExtinction, DataProducts.Extinction, DataProducts.ABackscatterCoeff] = ...
% % % %     FindHSRLParameters(DataProducts.BetaMProfile,                         ...
% % % %     Counts.CountRate{Map.Combined,1},                  ...
% % % %     Counts.CountRate{Map.Molecular,1},                 ...
% % % %     ones(1,size(Counts.CountRate{Map.Combined,1},2)),  ...
% % % %     ones(1,size(Counts.CountRate{Map.Molecular,1},2)), ...
% % % %     Altitude.RangeOriginal,                            ...
% % % %     JSondeData.ReceiverScaleFactor);
% % % % % Applying some data mask for the HSRL data
% % % % if Options.flag.mask_data == 1
% % % %     Combined_masked = Counts.CountRate{Map.Combined,1};
% % % %     Combined_masked(Counts.CountRate{Map.Combined,1} < 1/(R_size/PulseInfo.BinWidth)) = nan;
% % % %     Combined_masked(Counts.CountRate{Map.Molecular,1} < 1/(R_size/PulseInfo.BinWidth)) = nan;
% % % %     DataProducts.ABackscatterCoeff(isnan(Combined_masked))          = nan;
% % % %     DataProducts.AExtinction(isnan(DataProducts.ABackscatterCoeff)) = nan;
% % % %     clear Combined_masked
% % % % end
% % % % clear R_size


end