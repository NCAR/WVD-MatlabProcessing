% Written By: Robert Stillwell
% Written On: October 16, 2020
% Written For: National Center for Atmospheric Research

function [Sp] = ReadPCASpectra(Paths,WaveL,Op)
%
% Inputs: Paths: A structure containing all of the file paths needed for
%                analysis
%         WaveL: A structure containing all the wavelengths of the lasers
%                used for analysis
%         Op:    A structure containing the full user-defined options
%
% Outputs: Sp:   A cell array containing all the information needed to
%                reconstruct spectra using principle component analysis
%
%% Loading data by looping over file base names and modifiers
for n=1:1:length(Paths.PCA.Spectra)
    % Checking for the center wavelengths of each trained file
    Files = FindPossibleFiles(Paths.PCASpec,Paths.PCA.Spectra{n});
    FileWL = FindCenterWaveLength(Files);
    for m=1:1:length(Paths.PCA.Wavelengths)
        % Checking how far the trained data is from the desired wavelength
        WaveDiffs = abs(FileWL - median(WaveL.(Paths.PCA.Wavelengths{m,1}).Value,'omitnan').*1e9);
        % Loading the nearest data set
        CWLogging(['    Loading: ',Files(WaveDiffs == min(WaveDiffs)).name,'\n'],Op,'Main')
        Sp.(Paths.PCA.Wavelengths{m,1}).(Paths.PCA.SpectraLabels{n,1}) =...
          load(fullfile(Paths.PCASpec,Files(WaveDiffs == min(WaveDiffs)).name));
    end
end
end

function [Files] = FindPossibleFiles(BasePath,Type)
%
% Inputs: BasePath: The full file of the PCA files to be loaded
%         Type:     The type of file ot load (letters before first period)
%
% Outputs: Files:   A list of possible PCA files matching this request
%
%%
Files = dir(fullfile(BasePath,[Type,'.*.mat']));
end

function [FileWL] = FindCenterWaveLength(Files)
%
% Inputs: Files:   A list of possible files to load
%
% Outputs: FileWL: The center wavelength from each file
%
%%
FileWL = str2num(cell2mat(arrayfun(@(A) A.name(end-11:end-4),Files,'Uni',false))); %#ok<ST2NM>
end