% Written by: Robert Stillwell 
% Written for: NCAR
% Modificication Info: Created August 1, 2018

function [Struct] = RecursiveCell2Struct(Cell,FieldNames)
%                     
% Inputs: Cell:       A cell array to be converted to a structure
%         FieldNames: The names of the variables in each cell
%                     
% Outputs: Struct:    A structure with the same information as the input
%                     cell array simply recast for readability
%                     
%% Recursively converting cell to strucutre
for m=1:1:size(Cell,1)
    if iscell(Cell{m,1})
        % Dive further into cell to convert sub-cells to sub-structures
        if size(FieldNames,2) == 2
            Cell{m,1} = RecursiveCell2Struct(Cell{m,1},FieldNames{m,2});
        end
    end
   FieldNames2{m,1} = FieldNames{m,1};  %#ok<AGROW>
end
Struct = cell2struct(Cell,FieldNames2);
end