


function [Sigma] = BuildHitranSpectrum(Lambda,Pressure,Temperature,Gas)
%
%
%
%
%% Making wavelength list
% Loading Hitran data
cd HitranData/; 
formatSpec = '%*s%*s%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fid = fopen('O2_HITRAN2012_760_781.txt');
Hitran = textscan(fid, formatSpec,'headerLines', 14,'Delimiter',',');
Hitran{end} = [];
Hitran = cell2mat(Hitran);
cd ..
% Converting wavelength to cm^-1
Nu = (1e7)./Lambda;

%%
Sigma = BuildAbsorptionSpectrum(Hitran,Nu,Pressure,Temperature,Gas);
end