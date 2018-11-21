% Written By: Robert Stillwell
% Written For: NCAR ASP
% Written On: November 19, 2018
% This function rebuilds spectra using a trained principle component
% analysis method. Principle components are calculated from a broad
% training set. Weights for those principle components are calculated as a
% function of temperature and pressure then surface polynomials are fit to
% the weight contours. Here rebuilding is done using those polynomial
% contours and the training set information. Note multiple temperatures are
% pressures can be fit at once but they must be passed as a column vector.

function [RebuildPoly] = Rebuild(MeanSpectrum,PolyFitParams,PrincipleComponents,PressRebuild,MeanP,SigmaP,TempRebuild,MeanT,SigmaT)
%
% Inputs: MeanSpectrum:        The mean spectrum of the entire training set
%                              of values at each wavelength
%         PolyFitParams:       An array of structures containing the
%                              coefficients needed by the polyvaln function
%                              to recreate the polynomial cotours
%         PrincipleComponents: An array of principle components from the
%                              training set at the same wavelengths and
%                              resolutions as the mean spectrum
%         PressRebuild:        A column array of pressures at which the
%                              spectrum of interest is to be rebuilt
%         MeanP:               The mean pressure of the training set
%         SigmaP:              The standard deviation of the pressures of
%                              the training set 
%         TempRebuild:         A column array of temperatures at which the
%                              spectrum of interest is to be rebuilt
%         MeanT:               The mean temperature of the training set
%         SigmaT:              The standard deviation of the temperatures
%                              of the training set 
%
% Outputs: RebuildPoly:        An array of rebuilt spectra. The rows are at
%                              the same resolution as the mean spectrum of
%                              the training set and the columns will
%                              correspond to the size of the temperature
%                              and pressure arrays
%
%% Constants (number of principle components used to rebuild spectrum)
PC2Use       = 25;
%% Function handles needed
NormalizePolyInputs   = @(T1,Tm,Sigma)    (T1-Tm)./Sigma;
%% Rebuilding spectrum from trained polynomials
PolyWeights2 = zeros(PC2Use,size(TempRebuild,1)); % Pre-allocating memory
for m=1:1:PC2Use
    % Rebuilding absorption spectrum weights with polynomial fits
    PolyWeights2(m,:) = polyvaln(PolyFitParams{m,1},[NormalizePolyInputs(TempRebuild,MeanT,SigmaT),NormalizePolyInputs(PressRebuild,MeanP,SigmaP)]);
end
% Multiplying principle components by the weights
RebuildPoly = PrincipleComponents(:,1:PC2Use)*PolyWeights2(1:PC2Use,:);
% Adding the mean spectrum back in to rebuild the full spectrum
RebuildPoly = RebuildPoly+ MeanSpectrum';
end

function ypred = polyvaln(polymodel,indepvar)
% polyvaln: evaluates a polynomial model as a function of its variables
% usage: ypred = polyvaln(polymodel,indepvar)
%
% arguments: (input)
%  indepvar - (n x p) array of independent variables as columns
%        n is the number of data points to evaluate
%        p is the dimension of the independent variable space
%
%        IF n == 1, then I will assume there is only a
%        single independent variable.
%
%  polymodel - A structure containing a regression model from polyfitn
%        polymodel.ModelTerms = list of terms in the model
%        polymodel.Coefficients = regression coefficients
%  
%        Note: A polymodel can be evaluated for any set of
%        values with the function polyvaln. However, if you
%        wish to manipulate the result symbolically using my
%        own sympoly tools, this structure should be converted
%        to a sympoly using the function polyn2sympoly.
%
% Arguments: (output)
%  ypred - nx1 vector of predictions through the model.
%
%
% See also: polyfitn, polyfit, polyval, polyn2sympoly, sympoly
%
% Author: John D'Errico
% Release: 1.0
% Release date: 2/19/06

% get the size of indepvar
[n,p] = size(indepvar);
if (n == 1) && (size(polymodel.ModelTerms,2)==1)
  indepvar = indepvar';
  [n,p] = size(indepvar);
elseif (size(polymodel.ModelTerms,2)~=p)
  error 'Size of indepvar array and this model are inconsistent.'
end

% Evaluate the model
nt = size(polymodel.ModelTerms,1);
ypred = zeros(n,1);
for i = 1:nt
  t = ones(n,1);
  for j = 1:p
    t = t.*indepvar(:,j).^polymodel.ModelTerms(i,j);
  end
  ypred = ypred + t*polymodel.Coefficients(i);
end
end

