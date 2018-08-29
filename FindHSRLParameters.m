

function [AExtinction, Extinction, ABackscatterCoeff] = FindHSRLParameters(BetaM,Combined,Molecular,OverlapC,OverlapM,Range,ReceiverScaleFactor)
%
%
%
%
%
%
%
%
%% Loading overlap corrections to overwrite inputs
OverlapM = ones(1,size(Combined,2));
OverlapC = ones(1,size(Combined,2));
 
%% Aerosol backscatter coefficient
% Finding the differential overlap correction between channels
DiffOverlap        = OverlapC./OverlapM;
% Finding the backscatter ratio
BackscatterRatio   = Combined./Molecular;
% Converting the molecular backscatter vector to a 2 dimensional array
BetaMProfile2D     = repmat(BetaM, size(BackscatterRatio,1), 1);
% Calculating the aerosol backscatter coefficient
ABackscatterCoeff  = ((BackscatterRatio.*ReceiverScaleFactor.*DiffOverlap)-1).*BetaMProfile2D;
%% Aerosol backscatter coefficient error
 
%% Total Extinction (calculated with molecular signal)
% Calculating hte range resolution 
DeltaR     = Range(2) - Range(1);
% Approximate total extinction of signal using the molecular channel
Extinction = (1./2./DeltaR).*log((OverlapM./circshift(OverlapM,1))    .* ...
                                 (BetaM./circshift(BetaM,1))          .* ...
                                 ((Range.^2)./((Range+DeltaR).^2))       .* ...
                                 (circshift(Molecular,[0,1])./Molecular));
% Getting rid of very low SNR data
Extinction(Molecular < 0.01) = nan;
%% Aerosol Extinction
% Approximate aerosol extinction 
AExtinction = Extinction - (8.*pi./3).*BetaM;
 
%% Total Extinction error
 
 
end