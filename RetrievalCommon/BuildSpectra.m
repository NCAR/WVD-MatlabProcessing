% Written By: Robert Stillwell
% Written On: December 2, 2020
% Written For: National Center for Atmospheric Research

function [S] = BuildSpectra(PCASpectra,T,P,Lambda,Options)
%
% Inputs: PCASpectra: Structure containing needed PCA strutures to rebuild
%                     any arbitrary spectrum. Top level should define the
%                     wavelength and the next level what spectra to rebuild
%                     per wavelength.
%         T:          2d array of temperatures to rebuild spectra at. Units
%                     should be in Kelvin. 
%         P:          2d array of pressures to rebuild spectra at. Units
%                     should be in atmospheres.
%         Lambda:     Structure containing wavelength data for desired
%                     rebuilt spectra
%         Options:    Structure containing user defined options
%                     
% Outputs: S:         A structure containing the same top and mid levels as 
%                     PCASpectra structure but with the bottom level being
%                     3d arrays containing the desired spectra. The array
%                     size is the same as T and P for the first 2 indecies
%                     with the 3rd index being wavelength. 
%
%% Determining what spectra to rebuild
[Cell,FN,~,~] = RecursiveStruct2Cell(PCASpectra);
%% Reshaping temperature and pressure arrays and storing in a cell array
B = cellfun(@(S) reshape(S.Value,length(S.Range)*length(S.TimeStamp),1),{T;P},'Uni',false);
%% Rebuilding spectra at all temperatures
for m=1:1:size(Cell)          % Looping over wavelengths
    for n=1:1:size(Cell{m,1}) % Looping over spectra types
        tic
        % Unpacking PCA data structure
        PCAData = PCASpectra.(FN{m,1}).(FN{m,2}{n,1});
        % Making PCA spectra
        Spectra = PCALineFitting(PCAData,B{2,1},B{1,1});
        % Reshaping spectra to a 3d array (alt,time,wavelength)
        S.(FN{m,1}).(FN{m,2}{n,1}) = reshape(Spectra',length(T.Range),length(T.TimeStamp),length(PCAData.Lambda));
        % Interpolating to transmission wavelength and setting units
        if strcmp(FN{m,2}{n,1},'Absorption')
            S.(FN{m,1}).(FN{m,2}{n,1}) = S.(FN{m,1}).(FN{m,2}{n,1})./1e4;    %[m^2]
            S.(FN{m,1}).AbsorptionObserved = Interp2Lambda(S.(FN{m,1}).(FN{m,2}{n,1}),PCAData.Lambda,Lambda.(FN{m,1}).Value);
        end
        if strcmp(FN{m,2}{n,1},'RayleighBr')
            % Calculating a "delta function" aerosol profile
            AerosolFieldSize = size(S.(FN{m,1}).(FN{m,2}{n,1}));
            S.(FN{m,1}).Aerosol = zeros(AerosolFieldSize);
            S.(FN{m,1}).Aerosol(:,:,ceil(AerosolFieldSize(3)/2)) = 1;
            % Normalizing the RB spectra to have an integral of 1
            S.(FN{m,1}).(FN{m,2}{n,1}) = S.(FN{m,1}).(FN{m,2}{n,1})./trapz(S.(FN{m,1}).(FN{m,2}{n,1}),3);
        end
        CWLogging(['  Building: ',FN{m,1},' ',FN{m,2}{n,1},' took: ',num2str(toc),' [sec] \n'],Options,'Retrievals')
    end
end
end

function [Interp] = Interp2Lambda(Spectra,PCAL,Lam)
%
% Inputs: Spectra: A 3d array of possible spectra with the dimensions of:
%                  (altitude, time, wavelength)
%         PCAL:    The wavelength of the PCA reconstruction  (units are nm)
%         Lam:     The observed MPD wavelengths in time       (units are m)
%
% Outputs: Interp: Interpolated absorption coefficients to the obsereved
%                  wavelength with the dimensions of: (altitude, time)
%
%% Pre-allocating desired wavelength array
Interp = zeros(size(squeeze(Spectra(:,:,1))));
%% Looping over all elements of the array and interpolating
for a=1:1:size(Interp,1)      % Looping over range
    % Interpolating in time to the observed wavelength
    Interp(a,:) = interp2((1:size(Spectra,2))',PCAL,squeeze(Spectra(a,:,:))',...
                          (1:size(Spectra,2))',Lam.*1e9);
end
end

function [CrossSection] = PCALineFitting(TrainingSet,PressRebuild,TempRebuild)
%
% Inputs: TrainingSet:        A structure containing all of the information
%                             needed to rebuild the water vapor absorption
%                             spectrum from the training set
%         PressRebuild:       A column array of pressures at which the
%                             spectrum of interest is to be rebuilt        
%         TempRebuild:        A column array of temperatures at which the
%                             spectrum of interest is to be rebuilt
%
% Outputs: CrossSection:      Units [cm^2]
%
%% Rebuilding full spectrum for each desired temperature and pressure
CrossSection = Rebuild(TrainingSet.MeanSpectrum, ...
                       TrainingSet.PolyFitParams,...
                       TrainingSet.PrincipleComponents, ...
                       PressRebuild ,TrainingSet.MeanP ,TrainingSet.SigmaP , ...
                       TempRebuild  ,TrainingSet.MeanT ,TrainingSet.SigmaT);
end

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
PC2Use       = 45;
%% Function handles needed
NormalizePolyInputs   = @(T1,Tm,Sigma) (T1-Tm)./Sigma;
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

% polyvaln: evaluates a polynomial model as a function of its variables
function ypred = polyvaln(polymodel,indepvar)
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
%% get the size of indepvar
[n,p] = size(indepvar);
if (n == 1) && (size(polymodel.ModelTerms,2)==1)
  indepvar = indepvar';
  [n,p] = size(indepvar);
elseif (size(polymodel.ModelTerms,2)~=p)
  error 'Size of indepvar array and this model are inconsistent.'
end
%% Evaluate the model
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%