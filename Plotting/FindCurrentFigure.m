% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August 16th, 2020

% If the current figure input is not set, this function will assume there
% are no other figures and number these figures accordingly. If there are
% other figures, this function returns the number of the maximum number
function [FigNum] = FindCurrentFigure
%
% Inputs: none
%
% Outputs: FigNum: The highest figure number currently availible
%%
Figures = double(findall(0,'type','figure'));
if isempty(Figures)
    FigNum = 0;
else
    FigNum = max(Figures);
end
end