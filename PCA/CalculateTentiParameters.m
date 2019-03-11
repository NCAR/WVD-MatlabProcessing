


function [X,Y,K,Nu0] = CalculateTentiParameters(Pressure,Temperature,FrequencyChange,Lambda,Const)
%
%
%
%
%
%
%% Calculating magnitude of the wave vector and thermal velocity 
K   = 4.*pi./Lambda./sin(Const.Theta/2);
Nu0 = sqrt(Const.Kb.*Temperature./Const.MAir);

%% Calculating range of X and Y for training
X = 2.*pi.*FrequencyChange./sqrt(2)./K./Nu0;
Y = Pressure./sqrt(2)./K./Nu0./Const.Viscosity;

end