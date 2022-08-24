% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August, 2022

function [Const,Counts,Data1D,Scan,Spectra,Surface,Possible] = LoadAndPrepDataForRetrievals(As,Chan,Cal,Data,Op,Options,Paths)
%
%
%
%
%
%
%% Loading data needed for processing
Const = DefineConstants;
%% Parsing out data from complete data files to more useful form
[Counts,Data1D,Scan,Possible] = IdentifyNeededInfo(Data,Cal,As,Chan,Const);
if not(Possible) % If all data doesn't exist, return empty data
   Counts = []; Data1D = [];Scan = [];Spectra = []; Surface = [];
   return 
end
%% Getting approximate atmosphere info or using default
try   % Try using onboard weather station to approximate atmosphere
   Surface             = Data.TimeSeries.WeatherStation;
   Surface.TimeStamp   = Surface.TimeStamp.*60.*60;
catch % Using a default atmosphere if on can not be found
   Surface.TimeStamp   = Op.WV.TimeStamp;
   Surface.Temperature = ones(size(Surface.TimeStamp)).*15; % [Celcius]
   Surface.Pressure    = ones(size(Surface.TimeStamp)).*841;% [Millibar]
end
%% Downsample and interpolate ancillary data to known MPD grid
Data1D  = RecursivelyInterpolate1DStructure(Data1D,Options.TimeStamp,'linear');
Surface = RecursivelyInterpolate1DStructure(Surface,Options.TimeStamp,'linear');
%% Converting weather station data to the correct units
Surface.Temperature = Surface.Temperature + Const.C2K;  % [Kelvin]
Surface.Pressure    = Surface.Pressure./Const.MBar2Atm; % [Atmospheres]
%% Loading data needed for processing
Spectra.PCA = ReadPCASpectra(Paths,Data1D.Wavelength,Op);
end
