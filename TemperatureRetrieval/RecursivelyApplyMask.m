% Written by: Robert Stillwell 
% Written for: NCAR
% Modificication Info: Created February 13, 2018

function [ReturnDataStructure] = RecursivelyApplyMask(OriginalDataStructure)
%
% Inputs: OriginalDataStructure: Data structure to to interpolate
%         NewTime:               An array of monotonically increasing time
%                                stamps to interpolate data to 
%         Method:                Mthod of interpolation to be used
%
% Outputs: ReturnDataStructure:  A structure containing all data as before
%                                but with data interpolated to desired grid
%
%% Converting the surface weather structure into a cell array
[Cell,FieldNames,~] = RecursiveStruct2Cell(OriginalDataStructure);
%% Recursively performing an interpolation of the cell contents
CellDataNew = ApplyMask(Cell,FieldNames);
%% Convert the surface weather cell array back to a structure
ReturnDataStructure = RecursiveCell2Struct(CellDataNew, FieldNames);
end

% This function performs interpolation on all elements in a cell array
function [CellDataNew] = ApplyMask(CellData,FieldNames) 
%                        
% Inputs: CellData:      Cell array to recursively interpolate
%         FieldNames:    The names of the variables in each cell 
%                        
% Outputs: CellDataNew:  Cell array containing all data as before but with 
%                        processing masks applied
%                        
%% Recursively interpolating data contained within the cell array
Counter = 1;
for m=1:1:size(CellData,1)
   if iscell(CellData{m,1}) 
%        fprintf([FieldNames{m,1},'--------------------------\n'])
       % Need to dive down further into the cell array
       CellDataNew{m,1} = ApplyMask(CellData{m,1},FieldNames{m,2});                       %#ok<*AGROW>
   else
%        fprintf(['     ',FieldNames{m,1},'\n'])
       % The bottom of the tree so interpolate if it's not a single number
       if any(contains(FieldNames,'Mask'))
           Mask = CellData{contains(FieldNames,'Mask')==1,1};
%            if strcmp(FieldNames{m,1},'Mask') ~= 1
               if all(size(CellData{m,1}) == size(Mask))
                  CellData{m,1}(Mask==1) = nan; 
               end
               CellDataNew{Counter,1} = CellData{m,1};
               Counter = Counter + 1;
%            end
       else
          CellDataNew{m,1} = CellData{m,1}; 
       end
   end
end
end
