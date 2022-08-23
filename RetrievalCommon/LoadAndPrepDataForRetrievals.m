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
[Counts,Data1D,Scan,Possible] = IdentifyNeededInfo(Data,Cal,As,Chan);
if not(Possible)
   Counts = []; Data1D = [];Scan = [];Spectra = []; Surface = [];
   return 
end
%% Getting approximate atmosphere info or using default
try
   Surface             = Data.TimeSeries.WeatherStation;
   Surface.TimeStamp   = Surface.TimeStamp.*60.*60;
catch
   Surface.TimeStamp   = Op.WV.TimeStamp;
   Surface.Temperature = ones(size(Surface.TimeStamp)).*15;
   Surface.Pressure    = ones(size(Surface.TimeStamp)).*0.83;
end
%% Downsample and interpolate ancillary data to known MPD grid
Data1D  = RecursivelyInterpolate1DStructure(Data1D,Options.TimeStamp,'linear');
Surface = RecursivelyInterpolate1DStructure(Surface,Options.TimeStamp,'linear');
Surface.Temperature = Surface.Temperature + Const.C2K;
Surface.Pressure    = Surface.Pressure./Const.MBar2Atm;
%% Loading data needed for processing
Spectra.PCA = ReadPCASpectra(Paths,Data1D.Wavelength,Op);
end
