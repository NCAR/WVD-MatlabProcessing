% Written by: Robert Stillwell 
% Written for: NCAR
% Modificication Info: Created August 15, 2018

% This function compiles all the strucutre information into a cell array
% with all data pulled from arbitrary depth to the top level cells
function [AllData,TopNames] = FindDataStructures(Data)
%
% Inputs: Data:      A general structure containing lidar data with
%                    arbitrary depths
%
% Outputs: AllData:  A cell array containing all known data structures at
%                    the top level 
%          TopNames: A cell array containing the names of each data cell
%
%% Converting data structure to a cell
RawData   = struct2cell(Data);
FieldName = fieldnames(Data);
%% Diving into structure to find lowest levels & pulling data to the top
AllData  = {}; TopNames = {}; Counter = 1;  %Pre-allocating cell & counter
for m=1:1:size(RawData,1)
    % Finding data using a recursive search
    [TempCell,TempName] = RecursivelyFindData(RawData{m,1},FieldName{m,1});
    % Saving data into a single layer data structure 
    if isstruct(TempCell)
        AllData{Counter,1} = TempCell; %#ok<*AGROW>
        TopNames{Counter,1} = TempName;
        Counter = Counter + 1;
    else
        for n=1:1:size(TempCell,1)
            AllData{Counter,1} = TempCell{n,1};
            TopNames{Counter,1} = TempName{n,1};
            Counter = Counter + 1;
        end
    end 
end
end

% This function finds data contained within a single cell
function [Struct2Plot,FieldName] = RecursivelyFindData(Cell,FieldName)
%
% Inputs: Cell:         Contents of a section of cell array containing data
%         FieldName:    The field name generated from conversion of
%                       structure data to cell data
%
% Outputs: Struct2Plot: Cell array containing all data from the input
%                       cell brought to the top level
%          FieldName:   Cell array with all fieldnames of the cell data
%
%% Checking if structure has a timestamp field or if we need to dive deeper
if isfield(Cell,'TimeStamp')
    Struct2Plot = Cell;
else
    clear FieldName
    SubStruct = struct2cell(Cell);
    SubField  = fieldnames(Cell);
    for m=1:1:size(SubStruct,1)
        [Struct2Plot{m,1},FieldName{m,1}] = RecursivelyFindData(SubStruct{m,1},SubField{m,1}); 
    end
end
end