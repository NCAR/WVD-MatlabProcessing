% Written by: Robert Stillwell 
% Written for: NCAR
% Modificication Info: Created August 1, 2018

function [Cell,FieldNames,TimeStamps,Range] = RecursiveStruct2Cell(Struct)
%                      
% Inputs: Struct:      A structure to be converted to a cell array
%                      
% Outputs: Cell:       A cell array with all data as the input structure
%          FieldNames: The names of the variables in each cell 
%          TimeStamps: Time stamps at the bottom level of each structure
%                      
%% Converting structure to a cell array
Cell           = struct2cell(Struct);
FieldNames     = fieldnames(Struct);
[~,TimeStamps] = RecursivelyCheckIsField(Struct,'TimeStamp');
[~,Range]      = RecursivelyCheckIsField(Struct,'Range');
%% Recursively checking if sub-cells should also be converted
for m=1:1:size(Cell,1)
    if isstruct(Cell{m,1})
        % Dive down into structure to convert sub-structures to sub-cells
        [Cell{m,1},FieldNames{m,2},TimeStamps{m,1},Range{m,1}] = RecursiveStruct2Cell(Cell{m,1});
    end
end
end