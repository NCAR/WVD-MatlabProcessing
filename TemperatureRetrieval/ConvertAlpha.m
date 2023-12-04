



function [TCurrent,DTAll] = ConvertAlpha(Alpha,Const,Data1D,Data2D,Options,Surf,Spectra,Op,GuessLapse,StartCond,Paths)
%
%
%
%
%% Constants used to run this function
Tolerance  = 0.05;      % Threshold of average temperature change to exit
MaxDt      = 2;         % Threshold of temperature change per iteration

%% Pulling out the water vapor contour to modify the number density of O2
[IsField,Humidity] = RecursivelyCheckIsField(Data2D, {'MPD','Humidity'});
[IsField2,Atmo   ] = RecursivelyCheckIsField(Data2D, {'NCIP'});
if IsField && IsField2
    % Applying the humidity mask
    Humidity.Value(ceil(Humidity.Mask)==1 | Humidity.Value<0) = 0;
else
    [~,Humidity] = RecursivelyCheckIsField(Data2D, {'Onboard','WV'});
    [~,Atmo    ] = RecursivelyCheckIsField(Data2D, {'Guess'});
end

%% Preparing atmospheric parameters
% Converting from g/m^3 to molecules per m^3
Humidity.Value = Humidity.Value./Const.MWV/1000;
% Calculating the number of molecules period
A = Atmo.Pressure.Value.*Const.Atm2Pa.*Const.Av./Const.R./Atmo.Temperature.Value;
% Calculating the mixing ratio of water vapor
Qwv = Humidity.Value./A;


Humidity.Value(Qwv>=0.1) = nan;
Qwv(Qwv>=0.1) = nan;

%% Loading Hitran information if needed
if strcmp(Options.Method,'LBL')
    Hitran = LoadHitranData(Paths.HitranSpec);
end

%% Calculating constants
Gamma = Const.G0*Const.MolMAir/Const.R;
%% Iterating
% Seeding the loop with initial guesses of temperature and lapse rate
ConstProfile = (Options.Range').*(GuessLapse);
% Creating current temperature structure
TCurrent.TimeStamp = Atmo.Temperature.TimeStamp;
TCurrent.Range     = Atmo.Temperature.Range;
if strcmp(StartCond,'Cold')
    TCurrent.Value     = repmat(ConstProfile,1,size(Alpha,2))*0 + 240;
elseif strcmp(StartCond,'Warm')
    TCurrent.Value     = repmat(ConstProfile,1,size(Alpha,2))*0 + 330;
elseif strcmp(StartCond,'Bootstrap')
    TCurrent.Value     = Atmo.Temperature.Value; % Reset previously
else
    TCurrent.Value     = repmat(ConstProfile,1,size(Alpha,2))+Surf.Temperature.Value';
end
% Looping
DTAll = zeros(1,Options.TempIter).*nan;
for m=1:1:Options.TempIter
    % Calculating the change in temperature by either method
    if strcmp(Options.Method,'PCA')
        DeltaT = ConvertPCA(Alpha,Const,Spectra.PCA,Atmo.Pressure,...
                            TCurrent,Humidity.Value,Data1D.Wavelength,Op);
    else
        DeltaT = ConvertLineByLine(Alpha,Const,Spectra.PCA,Atmo.Pressure,...
                                   TCurrent,Qwv,Data1D.Wavelength,Op,GuessLapse,Surf,Gamma,Hitran);
    end
    % Limiting the gradient possible
    DeltaT(abs(DeltaT)>MaxDt) = sign(DeltaT(abs(DeltaT)>MaxDt)).*MaxDt;
    % Outputting temperature state
    TempDiffAvg = mean(mean(abs(DeltaT),'omitnan'),'omitnan');
    DTAll(m) = TempDiffAvg; 
    CWLogging(sprintf('      Mean dT: %4.3f\n',TempDiffAvg),Op,'Retrievals')
    % Updating the current temperature
    TCurrent.Value = TCurrent.Value + DeltaT;
    if abs(TempDiffAvg) <= Tolerance
        break
    end
%     % Plotting just to see whats going on
%     if mod(m,5)==1
%         figure(101);
%         subplot(ceil(Options.TempIter/5),1,floor(m/5)+1)
%         pcolor(TCurrent.TimeStamp./60./60,TCurrent.Range./1e3,DeltaT);
%         shading flat; colorbar; caxis([-MaxDt,MaxDt]); colormap(gca,redblue(64));
%         ylabel('Altitude [km]')
%     end
end
%% Removing data where the Dt is clearly not reached a stable point
TCurrent.Value(abs(DeltaT)==MaxDt) = nan;
end

function [HitranData] = LoadHitranData(FileName)
%
% Input: FileName
%
% Output: HitranData:
%
%% Loading Hitran data and converting
fid = fopen(FileName);
Data = textscan(fid,'%s','delimiter','\n');
fclose(fid);
Data = char(Data{:});
% Converting data
HitranData = zeros(size(Data,1),8);
for p=1:1:size(Data,1)
    % Reading hitran data (file needs padding with leading blank space)
    HitranStruct = ReadHitranData([' ',Data(p,:)]);
    % Reformatting raw data to expected format
    HitranData(p,:) = [HitranStruct.Nu      ,HitranStruct.S         ,HitranStruct.EinsteinA,...
                       HitranStruct.ABroadHW,HitranStruct.SBroadHW  ,HitranStruct.LowEnergy,...
                       HitranStruct.TempDep ,HitranStruct.PressShift];
end
end

% Data defined on: https://www.cfa.harvard.edu/hitran/formats.html
function [HT] = ReadHitranData(String)
%
% Inputs:
%
% Outputs:
%
%% Reading double inputs
HT.Nu         = String(4:15);
HT.S          = String(16:25);
HT.EinsteinA  = String(26:35);
HT.ABroadHW   = String(36:40);
HT.SBroadHW   = String(41:45);
HT.LowEnergy  = String(46:55);
HT.TempDep    = String(56:59);
HT.PressShift = String(60:67);
%% Converting all strings to doubles
FieldNames = fieldnames(HT);
HT = struct2cell(HT);
for m=1:1:length(HT)
   HT{m} = str2double(HT{m}); 
end
HT = cell2struct(HT,FieldNames);
end

function [DeltaT] = ConvertLineByLine(Alpha,Const,PCA,Pressure,Temperature,WV,Wavelength,Op,GuessLapse,Surf,Gamma,Hitran)
%
%
%
%
%% Calculating the weight of absorption of each line vs. line center
Nu = 1e-2./Wavelength.O2Online.Value';
for m=1:1:size(Hitran,1)
    Sigma(:,:,m) = BuildNormHitranSpectrum(Hitran(m,:),Nu,Pressure.Value,Temperature.Value);
end
%%
% SpecNew = BuildSpectra(PCA,Temperature,Pressure,Wavelength,Op);
LS = CalculateLineStrength(Const, Const.Eo./Const.H./Const.C./100, Const.O2LineS.*100, 296, Temperature.Value);
% LineShape = SpecNew.O2Online.AbsorptionObserved./LS;
LineShape = Sigma./LS;
% Calculating lapse rate of the current temperature profile
Lapse = ones(1,length(Surf.Temperature.Value)).*GuessLapse;
% Calculating simplifying constants
[C1,C2,C3] = CalculateConstants(Const,Surf,Gamma,Lapse,Temperature.Value);
% Calculating the update temperature
DeltaT = (Alpha - C1(:,:,1).*C2(:,:,1).*LineShape(:,:,1).*Const.QO2.*(1-WV) - ...
                  C1(:,:,2).*C2(:,:,2).*LineShape(:,:,2).*Const.QO2.*(1-WV))./...
                 (C1(:,:,1).*C2(:,:,1).*C3(:,:,1).*LineShape(:,:,1).*Const.QO2.*(1-WV) + ...
                  C1(:,:,2).*C2(:,:,2).*C3(:,:,2).*LineShape(:,:,2).*Const.QO2.*(1-WV));
% DeltaT = (Alpha - C1(:,:,1).*C2(:,:,1).*LineShape(:,:,1).*Const.QO2.*(1-WV))./...
%                  (C1(:,:,1).*C2(:,:,1).*C3(:,:,1).*LineShape(:,:,1).*Const.QO2.*(1-WV));
end

function [DeltaT] = ConvertPCA(Alpha,Const,PCA,Pressure,Temperature,WV,Wavelength,Op)
%
%
%
%
%%
% Number of total molecules per m^3 (based on pressure)
n = Pressure.Value.*Const.Atm2Pa./Const.Kb./Temperature.Value;
% Number of O2 molecules is (total - water) * O2 dry mixing ratio
n = (n - WV).*Const.QO2;
% Converting absorption coefficient to absorption cross section
SigmaM = Alpha./n;
% Calculating the absorption cross section at the current temperature
SpecNew = BuildSpectra(PCA,Temperature,Pressure,Wavelength,Op);
% Calculating the derivative of absorption cross section at the current temperature
SpecNewDer = BuildSpectraDerivative(PCA,Temperature,Pressure,Wavelength,Op);
% Calculating the change in temperature required
DeltaT = (SigmaM - SpecNew.O2Online.AbsorptionObserved)./ ...
                                  (SpecNewDer.O2Online.AbsorptionObserved);
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
LS = (LS0.*((RefT./Temperature).^1.0).*exp((100.*Const.H.*Const.C./Const.Kb) ...
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


