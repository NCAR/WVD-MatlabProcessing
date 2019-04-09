% Written by: Robert Stillwell
% Written for: National Center For Atmospheric Research
% This function processes all the High Spectral Resolution lidar retrievals 
% from the MPD system. It also filters and averages the data.
% Modification info: Created: November 13, 2018

function [] = RetrievalsHSRL(Altitude, Capabilities, Counts, DataProducts, Map, Options, Paths, PulseInfo, SurfaceWeather)
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

%% Spectrum parameters to rebuild
RebuildFreq  = linspace(-5e9,5e9,250);
CenterLam    = 770.1085e-9;     % Center wavelength to calculate/plot [nm]

RebuildPress = reshape(DataProducts.PresssureProfile,prod(size(DataProducts.PresssureProfile)),1); %#ok<*PSIZE>
RebuildTemp  = reshape(DataProducts.TemperatureProfile,prod(size(DataProducts.TemperatureProfile)),1);

% RebuildPolyQ = RebuildRBSpectra(RebuildPress,RebuildTemp,RebuildFreq,CenterLam,Paths);



% % % %% Calculate HSRL parameters
% % % R_size                = 150;
% % % % Calculating the HSRL retrievals
% % % [DataProducts.AExtinction, DataProducts.Extinction, DataProducts.ABackscatterCoeff] = ...
% % %     FindHSRLParameters(DataProducts.BetaMProfile,                         ...
% % %                        Counts.CountRate{Map.Combined,1},                  ...
% % %                        Counts.CountRate{Map.Molecular,1},                 ...
% % %                        ones(1,size(Counts.CountRate{Map.Combined,1},2)),  ...
% % %                        ones(1,size(Counts.CountRate{Map.Molecular,1},2)), ...
% % %                        Altitude.RangeOriginal,                            ...
% % %                        JSondeData.ReceiverScaleFactor);
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

function [RebuildPolyQ] = RebuildRBSpectra(RebuildPress,RebuildTemp,RebuildFreq,CenterLam,Paths)
%
%
%
%
%
%% Loading the training set
load(Paths.PCARBSet); 
%% Calculating rebuilt Tenti parameters
fprintf('      HSRL Retrieval: Calculating Rayleigh Brillouin Spectra\n')
[X,Y,~,~]   = CalculateTentiParametersNDim(RebuildPress,RebuildTemp,RebuildFreq,CenterLam,Const);
X           = shiftdim(X,size(size(X),2)-1);
RebuildPoly = RebuildNDim(MeanSpectrum,PolyFitParams,PrincipleComponents,Y);
%% Interpolating spectra back to desired frequency grid
RebuildPolyQ = zeros(size(X)); % Pre-allocating data storage array
for m=1:1:size(X,3)       % Looping over the columns
    for n=1:1:size(X,2)   % Looping over the rows
        % Interpolating retrieved spectrum to frequency grid of interest
        RebuildPolyQ(:,n,m) = interp1(XLimits',squeeze(RebuildPoly(:,n,m)),squeeze(X(:,n,m)));      
    end
end
% Normalizing all spectra
RebuildPolyQ = RebuildPolyQ./trapz(RebuildFreq./1e9,RebuildPolyQ);
end

function [X,Y,K,Nu0] = CalculateTentiParametersNDim(Pressure,Temperature,FrequencyChange,Lambda,Const)
%
%
%
%
%
%
%% Calculating magnitude of the wave vector and thermal velocity 
K   = 4.*pi./Lambda./sin(Const.Theta/2);
Nu0 = sqrt(Const.Kb.*Temperature./Const.MAir);
%% Reshaping frequency array to build on size of temperature array
FrequencyChange = shiftdim(repmat(FrequencyChange',fliplr([fliplr(size(Nu0)),1])),1);
%% Calculating range of X and Y for training
X = 2.*pi.*FrequencyChange./sqrt(2)./K./Nu0;
Y = Pressure./sqrt(2)./K./Nu0./Const.Viscosity;
end

function [RebuildPoly] = RebuildNDim(MeanSpectrum,PolyFitParams,PrincipleComponents,Y)
%
% Inputs: MeanSpectrum:        The mean spectrum of the entire training set
%                              of values at each wavelength
%         PolyFitParams:       An array of structures containing the
%                              coefficients needed by the polyvaln function
%                              to recreate the polynomial cotours
%         PrincipleComponents: An array of principle components from the
%                              training set at the same wavelengths and
%                              resolutions as the mean spectrum
%         Y:                   Array of Y parameter values to be rebuilt 
%
% Outputs: RebuildPoly:        An array of rebuilt spectra. The rows are at
%                              the same resolution as the mean spectrum of
%                              the training set and the columns will
%                              correspond to the size of the temperature
%                              and pressure arrays
%
%% Constants (number of principle components used to rebuild spectrum)
PC2Use       = 10;
%% Rebuilding spectrum from trained polynomials
RebuildPoly = zeros([size(PrincipleComponents,1),size(Y)]);
for m=1:1:PC2Use
    % Rebuilding absorption spectrum weights with polynomial fits
    PolyWeights2 = polyval(PolyFitParams{m,1},Y);
    % Multiplying by the weights by the principle components
    RebuildPoly = RebuildPoly + PrincipleComponents(:,m).*shiftdim(PolyWeights2,-1);
end
% Adding the mean spectrum back in to rebuild the full spectrum
RebuildPoly = RebuildPoly+ MeanSpectrum';
end






