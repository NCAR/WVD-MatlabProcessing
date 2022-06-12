



function [TCurrent,DTAll] = ConvertAlpha2Temperature(Alpha,Const,Data1D,Data2D,Options,Surf,Spectra,Op,GuessLapse,StartCond)
%
%
% Add moisture mixing ratio
% Add ability to check lapse rate
%
%
%% Constants used to run this function
Tolerance  = 0.05;      % Threshold of average temperature change to exit
MaxDt      = 2;         % Threshold of temperature change per iteration

Qwv = 0;

%% Calculating constants
Gamma = Const.G0*Const.MolMAir/Const.R;
%% Iterating
% Seeding the loop with initial guesses of temperature and lapse rate
ConstProfile = (Options.Range').*(GuessLapse);
% Creating current temperature structure
TCurrent.TimeStamp = Data2D.NCIP.Temperature.TimeStamp;
TCurrent.Range     = Data2D.NCIP.Temperature.Range;
if strcmp(StartCond,'Cold')
    TCurrent.Value     = repmat(ConstProfile,1,size(Alpha,2))*0 + 240;
elseif strcmp(StartCond,'Warm')
    TCurrent.Value     = repmat(ConstProfile,1,size(Alpha,2))*0 + 330;
elseif strcmp(StartCond,'Bootstrap')
    TCurrent.Value     = Data2D.NCIP.Temperature.Value; % Reset previously
else
    TCurrent.Value     = repmat(ConstProfile,1,size(Alpha,2))+Surf.Temperature.Value';
end
% Looping
for m=1:1:Options.TempIter
    % Calculating the absorption lineshape function (update with temp)
    PCA.O2Online.Absorption = Spectra.PCA.O2Online.Absorption;
    SpecNew = BuildSpectra(PCA,TCurrent,Data2D.NCIP.Pressure,Data1D.Wavelength,Op);
    LS = CalculateLineStrength(Const, Const.Eo./Const.H./Const.C./100, Const.O2LineS.*100, 296, TCurrent.Value);
    LineShape = SpecNew.O2Online.AbsorptionObserved./LS;
    % Calculating lapse rate of the current temperature profile 
%     Lapse = FittingLapseRate(TCurrent);
    Lapse = ones(1,length(Surf.Temperature.Value)).*GuessLapse; 
    % Calculating simplifying constants
    [C1,C2,C3] = CalculateConstants(Const,Surf,Gamma,Lapse,TCurrent.Value);
    % Calculating the update temperature
%     DeltaT = Alpha./(C1.*C2.*C3.*LineShape.*Const.QO2.*(1-Qwv)) - 1./C3;
    DeltaT = (Alpha - C1(:,:,1).*C2(:,:,1).*LineShape(:,:,1).*Const.QO2.*(1-Qwv) - ...
                      C1(:,:,2).*C2(:,:,2).*LineShape(:,:,2).*Const.QO2.*(1-Qwv))./...
                     (C1(:,:,1).*C2(:,:,1).*C3(:,:,1).*LineShape(:,:,1).*Const.QO2.*(1-Qwv) + ...
                      C1(:,:,2).*C2(:,:,2).*C3(:,:,2).*LineShape(:,:,2).*Const.QO2.*(1-Qwv));
    % Limiting the gradient possible
    DeltaT(abs(DeltaT)>MaxDt) = sign(DeltaT(abs(DeltaT)>MaxDt)).*MaxDt;
    % Outputting temperature state
    TempDiffAvg = mean(mean(abs(DeltaT),'omitnan'),'omitnan');
    DTAll(m) = TempDiffAvg; %#ok<AGROW>
    CWLogging(sprintf('      Mean dT: %4.3f\n',TempDiffAvg),Op,'Retrievals')
    % Updating the current temperature
    TCurrent.Value = TCurrent.Value + DeltaT;
    if abs(TempDiffAvg) <= Tolerance
        break
    end
% % %     % Plotting just to see whats going on
% % %     if mod(m,5)==1
% % %         figure(101);
% % %         subplot(ceil(Options.TempIter/5),1,floor(m/5)+1)
% % %         pcolor(TCurrent.TimeStamp./60./60,TCurrent.Range./1e3,DeltaT);
% % %         shading flat; colorbar; caxis([-MaxDt,MaxDt]); colormap(gca,redblue(64));
% % %         ylabel('Altitude [km]')
% % %     end
end
%% Removing data where the Dt is clearly not reached a stable point
TCurrent.Value(abs(DeltaT)==2) = nan;
end

function [LS] = CalculateLineStrength(Const, Energy, LS0, RefT, Temperature)
%
% Inputs: Const:       A structure containing needed universal constants in
%                      standard MKS units
%         Energy:      The energy of the selected transition with units of
%                      [cm^-1]
%         LS0:         The tabulated line strength of the transition from
%                      Hitran at standard temperature and pressure
%         RefT:        The reference pressure of the tabulated Lorentz
%                      broadening halfwidth with units of [K]
%         Temperature: The temperature of the desired Lorentz broadeding
%                      calculation with units of [K]
%
% Outputs: LS:         The line strength of the transition at the desired
%                      temperature (units of......)
%
%% Calculating the line strength
LS = (LS0.*((RefT./Temperature).^1.5).*exp((100.*Const.H.*Const.C./Const.Kb) ...
                                     .*Energy.*(1./RefT - 1./Temperature)))./100;
end


function [C1,C2,C3] = CalculateConstants(Const,Surf,Gamma,Lapse,Tc)
%
% Tc = temp current
%
%
%
%% Calculating intermediate variable for simplicity
GL = Gamma./Lapse;
%% Calculating the simplifying constant variables Kevin uses
% Updates each time because the lapse rate updates
C1 = Const.O2LineS.*Const.HitranTo.*(Surf.Pressure.Value'.*Const.Atm2Pa).* ...
     exp(Const.Eo/Const.Kb/Const.HitranTo)./(Const.Kb.*Surf.Temperature.Value'.^(-GL));
% Updates each time because it depends on current iteration's temperature 
C2 = Tc.^(-GL-2).*exp(-Const.Eo./Const.Kb./Tc);
% Updates each time because it depends on current iteration's temperature
C3 = (-GL - 2)./Tc + Const.Eo./Const.Kb./(Tc.^2);
end

function [LapseRates] = FittingLapseRate(TCurrent)
%
%
%
%% Setting fit perameters
RangeAllowed = [1e3,2.5e3];

Indices = [find(TCurrent.Range>RangeAllowed(1),1,'first'),...
           find(TCurrent.Range>RangeAllowed(2),1,'first')];

%% Pre-allocating data
LapseRates = zeros(1,length(TCurrent.TimeStamp));
%% Least squares fitting data to find the lapse rate
tic
for m=1:1:length(LapseRates)
    % Grabbing data to check for nans
    TData     = TCurrent.Value(Indices(1):Indices(2),m);
    RangeData = TCurrent.Range(Indices(1):Indices(2));
    % Removing nan values from data
    RangeData(isnan(TData)) = [];
    TData(isnan(TData)) = []; 
    % Making weight table (intercept, slope)
    X = [ones(size(TData)),RangeData'];
    % Doing a least squares fit
    Coeffs = X\TData;
    % Saving data
    LapseRates(m) = Coeffs(2);
end
fprintf(['      Lapse rate fitting took: ',num2str(toc),' [sec] with mean: ',num2str(mean(LapseRates)),'\n'])
end

