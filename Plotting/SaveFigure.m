% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created March 29th, 2021

function SaveFigure(FN,Op,Paths,Type)
%
% Inputs: FN:    The figure number to save and close
%         Op:    A structure containing all of the user defined options
%         Paths: A structure containing all the user defined file path info
%         Type:  A string containing the type of figure to be saved
%
%%
if Op.SaveFigures
    Base = fullfile(Paths.Quicklook,Type);
    if ~exist(Base,'dir')
        mkdir(Base)
    end
    saveas(FN,fullfile(Base,[lower(erase(Op.System,'_')),'.',Op.Date,'.',Type,'.png']))
    %saveas(FN,fullfile(Base,[lower(erase(Op.System,'_')),'.',Op.Date,'.',Type,'.fig']))
    close(FN)
end
end
