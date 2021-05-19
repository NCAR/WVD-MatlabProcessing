% Written By: Robert Stillwell
% Written On: October 16, 2020
% Written For: National Center for Atmospheric Research

function [Sp] = ReadPCASpectra(Paths)
%
% Inputs: Paths: A structure containing all of the file paths needed for
%                analysis
%
% Outputs: Sp:   A cell array containing all the information needed to
%                reconstruct spectra using principle component analysis
%
%% Moving to path containing PCA data
cd(Paths.PCASpec)  
%% Loading data by looping over file base names and modifiers
for m=1:1:length(Paths.PCA.Wavelengths)
    for n=1:1:length(Paths.PCA.Spectra)
        Sp.(Paths.PCA.Wavelengths{m,1}).(Paths.PCA.SpectraLabels{n,1}) =...
          load([Paths.PCA.Wavelengths{m,1},Paths.PCA.Spectra{n,1},'.mat']);
    end
end
%% Moving to path containg code
cd(Paths.Code)     
end