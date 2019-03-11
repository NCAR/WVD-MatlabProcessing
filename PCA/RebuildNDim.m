% Written By: Robert Stillwell
% Written For: NCAR ASP
% Written On: December 2, 2018
% This function rebuilds spectra using a trained principle component
% analysis method. Principle components are calculated from a broad
% training set. Weights for those principle components are calculated as a
% function of temperature and pressure then surface polynomials are fit to
% the weight contours. Here rebuilding is done using those polynomial
% contours and the training set information.

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
% % PolyWeights2 = shiftdim(zeros([PC2Use,size(Y)]),0); % Pre-allocating memory
RebuildPoly = zeros([size(PrincipleComponents,1),size(Y)]);

for m=1:1:PC2Use
    % Rebuilding absorption spectrum weights with polynomial fits
    PolyWeights2 = polyval(PolyFitParams{m,1},Y);
    % Multiplying by the weights by the principle comonents
    RebuildPoly = RebuildPoly + PrincipleComponents(:,m).*shiftdim(PolyWeights2,-1);
end
% % Multiplying principle components by the weights
% RebuildPoly = PrincipleComponents(:,1:PC2Use)*PolyWeights2(1:PC2Use,:);
% Adding the mean spectrum back in to rebuild the full spectrum
RebuildPoly = RebuildPoly+ MeanSpectrum';
end


