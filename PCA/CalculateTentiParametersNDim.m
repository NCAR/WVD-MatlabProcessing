




function [X,Y,K,Nu0] = CalculateTentiParametersNDim(Pressure,Temperature,FrequencyChange,Lambda,Const)
%
%
%
%
%
%
%% Calculating magnitude of the wave vector and thermal velocity 
K   = 4.*pi./Lambda./sin(Const.Theta/2);
Nu0 = sqrt(Const.Kb.*Temperature./Const.MAir);

%% Reshaping frequency array to build on size of temperature array
FrequencyChange = shiftdim(repmat(FrequencyChange',fliplr([fliplr(size(Nu0)),1])),1);

%% Calculating range of X and Y for training
X = 2.*pi.*FrequencyChange./sqrt(2)./K./Nu0;
Y = Pressure./sqrt(2)./K./Nu0./Const.Viscosity;

end