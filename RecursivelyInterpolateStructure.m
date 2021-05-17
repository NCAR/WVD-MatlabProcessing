% Written by: Robert Stillwell 
% Written for: NCAR
% Modificication Info: Created February 13, 2018

function [ReturnDataStructure] = RecursivelyInterpolateStructure(OriginalDataStructure,NewTime,Method)
%
% Inputs: OriginalDataStructure: Data structure to to interpolate
%         NewTime:               An array of monotonically increasing time
%                                stamps to interpolate data to 
%         Method:                Mthod of interpolation to be used
%
% Outputs: ReturnDataStructure:  A structure containing all data as before
%                                but with data interpolated to desired grid
%
%% Converting the structure into a cell array
[Cell,FieldNames,TimeStamps] = RecursiveStruct2Cell(OriginalDataStructure);
%% Recursively performing an interpolation of the cell contents
CellDataNew = RecursiveInterpolateData(Cell,TimeStamps,NewTime,Method,FieldNames);
%% Convert the cell array back to a structure
ReturnDataStructure = RecursiveCell2Struct(CellDataNew, FieldNames);
end

% This function performs interpolation on all elements in a cell array
function [CellDataNew] = RecursiveInterpolateData(CellData,OldTime,NewTime,Method,FieldNames) 
%                        
% Inputs: CellData:      Cell array to recursively interpolate
%         OldTime:       Array of time stamps to interpolate from
%         NewTime:       Array of time stamps to interpolate data to 
%         Method:        The method of interpolation. Inputs should be the
%                        same as those for the interp1 function
%                        
% Outputs: CellDataNew:  Cell array containing all data as before but with 
%                        data interpolated to desired grid
%                        
%% Recursively interpolating data contained within the cell array
for m=1:1:size(CellData,1)
   if iscell(CellData{m,1}) 
       % Need to dive down further into the cell array
% %        fprintf([FieldNames{m,1},'--------------------------\n'])
       Temp = RecursiveInterpolateData(CellData{m,1},OldTime{m,1},NewTime,Method,FieldNames{m,2});
       CellDataNew{m,1} = Temp;                       %#ok<*AGROW>
   else
       % The bottom of the tree so interpolate if it's not a single number
% %        fprintf(['     ',FieldNames{m,1},'\n'])
       if size(OldTime,1) == size(CellData{m,1},1)
           if length(OldTime) >2 % Checking to see if interpolation will fail
               if size(OldTime,1) ~= size(NewTime,1) % Already interpolated
                   CellData{m,1} = interp1(OldTime,CellData{m,1},NewTime,Method);
               end
           else
               % Filling data with default status (can't interpolate 0/1 data point)
               if strcmp(FieldNames{m,1},'TimeStamp'); CellData{m,1} = NewTime;
               else;                                   CellData{m,1} = NewTime.*nan;
               end
           end
       end
       CellDataNew{m,1} = CellData{m,1};
   end
end
end
