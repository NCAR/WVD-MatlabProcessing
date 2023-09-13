% Written By: Robert Stillwell
% Written On: December 4, 2020
% Written For: National Center for Atmospheric Research

function [Alpha,Pt] = PerturbativeRetrieval(Constants,Counts,Data2D,Options,Spectra,Op)
%
%
%
%
%
%
%% Checking which data inputs to use
if strcmp(Options.HSRLType,'Py') || strcmp(Options.HSRLType,'PyP')
    T    = Data2D.NCIP.Temperature.Value;
    P    = Data2D.NCIP.Pressure.Value;
    HSRL = Data2D.MPD.BackRatio.Value;
    WV   = Data2D.MPD.Humidity.Value;
else
    T    = Data2D.Guess.Temperature.Value;
    P    = Data2D.Guess.Pressure.Value;
    HSRL = Data2D.Onboard.HSRL.Smoothed;
    WV   = Data2D.Onboard.WV.Value;
end
%% Zeroth order retrieval
tic
Pt.Order0 = CalculateAlpha0(Counts.BGSub.O2Online.Counts,          ...
                            Counts.BGSub.O2Offline.Counts,         ...
                            Options.BinRange,                      ...
                            P,                               ...
                            T,                               ...
                            WV,Constants,   ...
                            Options,Spectra.Rebuilt);
CWLogging(['  Calculating: Alpha0 took: ',num2str(toc),' [sec] \n'],Op,'Retrievals')
%% First order retrieval
tic
Pt.Order1 = CalculateAlpha1(Pt.Order0,Spectra,HSRL,Options.BinRange,Options.Range');
CWLogging(['  Calculating: Alpha1 took: ',num2str(toc),' [sec] \n'],Op,'Retrievals')
%% Second order retrieval
tic
Pt.Order2 = CalculateAlpha2(Pt.Order0,Pt.Order1,Spectra);
CWLogging(['  Calculating: Alpha2 took: ',num2str(toc),' [sec] \n'],Op,'Retrievals')
Alpha = Pt.Order0.Alpha.O2Online + Pt.Order1.Alpha.O2Online + Pt.Order2.Alpha.O2Online;
end

function [Order0] = CalculateAlpha0(Online,Offline,DeltaR,P,T,WV,Const,Options,Spectra)
%
% Inputs: Online:  Background subtracted counts in the online channel
%         Offline: Background subtracted counts in the offline channel
%         DeltaR:  Range resolution of count data                   [m]
%         P:       Estimated pressure of the atmosphere             [atm] 
%         T:       Estimated temperature of the atmosphere          [K]
%         WV:      Observed water vapor number density              [g/m^3]
%         Const:   A structure containing all of the universal and earth
%                  constants needed for temperature analysis
%         Options:
%         Spectra: 
%
% Outputs: Order0: A structure containing the following variables:
%          -Alpha: Zeroth order absorption coefficient              [m^-1]
%          -Tm:    Spectral transmission of the light as a function of
%                  (alt,time,frequency)                             [unitless] 
%
%% Standard DIAL derivative term
A = (Online.*circshift(Offline,[1,0]))./(circshift(Online,[1,0]).*Offline);
A = circshift(A,[-1,0]); % Removing offset of 1 bin upwards shift 
A = A(1:size(P,1),:);
%% Number density of O2 (Total air - moisture = dry * QO2 = dry O2)
WV(WV<0) = 0;
NO2 = (P.*Const.Atm2Pa./Const.Kb./T - WV./(1e3.*Const.MWV)).*Const.QO2;
%% Alpha0 and transmission from alpha0
Order0.Alpha.O2Offline = Spectra.O2Offline.AbsorptionObserved.*NO2;
Order0.Alpha.O2Online  = real(Order0.Alpha.O2Offline - log(A)./2./DeltaR);
%% Looping over effective absorption to calculate further needed fields
[~,FN] = RecursiveStruct2Cell(Order0.Alpha);
for m=1:1:length(FN)
    % Calculating transmittion as function of (alt,time,frequency)...note
    % abs. array, F, is normalized because alpha used will belong to each
    F = Spectra.(FN{m}).Absorption;
    F = F./max(F,[],3);
    % Removing points where WV or A was not known
    NanMap = (isnan(WV) | isnan(Order0.Alpha.(FN{m})));
    Order0.Alpha.(FN{m})(NanMap) = 0;
    % Calculating spectral transmission
    Order0.Tm.(FN{m}) = exp(-cumtrapz(Options.Range,Order0.Alpha.(FN{m}).*F,1)); 
    % Removing points that can not be known because no WV info exists
    Order0.Alpha.(FN{m})(NanMap) = nan;
    Order0.Tm.(FN{m})(NanMap)    = nan;
end
end

function [Order1] = CalculateAlpha1(Order0,Spectra,BSR,DeltaR,Range)
%
% Inputs: Order0:
%         Spectra:
%         BSR:
%         DeltaR:
%
%% Filtering out bad data and calculating the inverse backscatter ratio
BSR(BSR<1) = 1;
IBSR = 1./BSR;
%% Calculating supporting variables
[~,FieldNames] = RecursiveStruct2Cell(Order0.Alpha);
for m=1:1:length(FieldNames)
    % Normalized absorption line shape (normalization needs to be to the  
    % online max so that spectra are on the same scale and normalized to 1)
    Order1.f.(FieldNames{m}) = Spectra.Rebuilt.(FieldNames{m}).Absorption./max(Spectra.Rebuilt.O2Online.Absorption,[],3);
    % Receiver optical transmission
    E = permute(Spectra.Optics.(FieldNames{m}).Transmission,[1,3,2]);
    % Lineshape of the backscattered light
    g = (1-IBSR).*Spectra.Rebuilt.(FieldNames{m}).Aerosol + ...
        IBSR .*Spectra.Rebuilt.(FieldNames{m}).RayleighBr;
    % Derivative in range of g    (dimensions are (alt, time, frequency))
    % ....might be negative of what I want
    dgdr = (g - circshift(g,[1,0,0]))./DeltaR;
    % Temperary variables to simplify the math
    Order1.Eta.(FieldNames{m})  = dgdr.*E.*Order0.Tm.(FieldNames{m});
    Order1.Zeta.(FieldNames{m}) = g.*E.*Order0.Tm.(FieldNames{m});
    % Calculate perturbative terms
    Order1.DeltaW.(FieldNames{m}) = trapz(Order1.Zeta.(FieldNames{m}).*(1-Order1.f.(FieldNames{m})),3)./trapz(Order1.Zeta.(FieldNames{m}),3);
    Order1.DeltaG.(FieldNames{m}) = trapz(Order1.Eta.(FieldNames{m}),3)./trapz(Order1.Zeta.(FieldNames{m}),3);
end
clear f E g dgdr
%% Calculating the 1st order absorption coefficient & transmission
for m=1:1:length(FieldNames)
    % Calculating alpha
    Order1.Alpha.(FieldNames{m}) = (Order0.Alpha.(FieldNames{m}).*Order1.DeltaW.(FieldNames{m}) + Order1.DeltaG.O2Online - Order1.DeltaG.O2Offline)./2;
    % Getting rid of nan values so intergration is possible
    IsNanMap = isnan(Order1.Alpha.(FieldNames{m}));
    Order1.Alpha.(FieldNames{m})(IsNanMap) = 0;
    % Calculating spectral transmission 
    Order1.Tm.(FieldNames{m}) = exp(-cumtrapz(Range,Order1.Alpha.(FieldNames{m}).*Order1.f.(FieldNames{m}),1));
    Order1.Tm.(FieldNames{m})(IsNanMap) = nan;
    % Reapplying nan values 
    Order1.Alpha.(FieldNames{m})(IsNanMap) = nan;
end
end

function [Order2] = CalculateAlpha2(Order0,Order1,Spectra)
%
%
%
%
%% Calculating 2nd order parameters for online and offline
[~,FN] = RecursiveStruct2Cell(Order0.Alpha);
for m=1:1:length(FN)
    % Normalized absorption line shape (normalization needs to be to the  
    % online max so that spectra are on the same scale and normalized to 1)
    f = Spectra.Rebuilt.(FN{m}).Absorption./max(Spectra.Rebuilt.O2Online.Absorption,[],3);
    % Calculating the integral of Zeta as it is used a lot
    IntZeta = trapz(Order1.Zeta.(FN{m}),3);    
    % Calculate perturbative terms 
    Order2.DeltaG.(FN{m}) = Order1.DeltaG.(FN{m}).*trapz(Order1.Zeta.(FN{m}).*(1-Order1.Tm.(FN{m})),3)./IntZeta - ...
                               trapz(Order1.Eta.(FN{m}).*(1-Order1.Tm.(FN{m})),3)./IntZeta; 
    Order2.DeltaW.(FN{m}) = trapz(Order1.Zeta.(FN{m}).*(1-f),3).*trapz(Order1.Zeta.(FN{m}).*(1-f).*(1-Order1.Tm.(FN{m})),3)./(IntZeta.^2) - ...
                               trapz(Order1.Zeta.(FN{m}).*(1-f).*(1-Order1.Tm.(FN{m})),3)./IntZeta;
    
end
%% Calculating the second order absoprtion 
for m=1:1:length(FN)
    Order2.Alpha.(FN{m}) = (Order1.Alpha.(FN{m}).*Order1.DeltaW.(FN{m}) + ...
                            Order0.Alpha.(FN{m}).*Order2.DeltaW.(FN{m}) + ...
                            Order2.DeltaG.O2Online - Order2.DeltaG.O2Offline)./2;
end
end