% Written by: Robert Stillwell 
% Written for: NCAR
% Modificication Info: Created February 13, 2018

function [ReturnDataStructure] = RecursivelyInterpolate1DStructure(OriginalDataStructure,NewTime,Method)
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
[Cell,FieldNames,TimeStamps] = RecursiveStruct2Cell(OriginalDataStructure);
%% Recursively performing an interpolation of the cell contents
CellDataNew = RecursiveInterpolateData(Cell,TimeStamps,NewTime,Method,FieldNames);
%% Convert the surface weather cell array back to a structure
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
%        fprintf([FieldNames{m,1},'--------------------------\n'])
       % Need to dive down further into the cell array
       Temp = RecursiveInterpolateData(CellData{m,1},OldTime{m,1},NewTime,Method,FieldNames{m,2});
       CellDataNew{m,1} = Temp;                       %#ok<*AGROW>
   else
       % The bottom of the tree so interpolate if it's not a single number
%        fprintf(['     ',FieldNames{m,1},'\n'])
       if all(size(OldTime(:)) == size(CellData{m,1}(:)))
           if size(OldTime,1) > 1 % Does more than 1 data point exists
               % Downsample data
               A = OldTime(:); B = NewTime(:);
               Downsample = ceil((B(2)-B(1))/(A(2)-A(1)));
               % Interpolating the downsampled data
               if sum( isnan(OldTime(:))) > 0 
                   CellData{m,1}(isnan(OldTime(:))) = [];
                   IntTime = OldTime;
                   IntTime(isnan(IntTime)) = [];
               else
                   IntTime = OldTime;
               end
               
               CellData{m,1} = interp1(movmean(IntTime(:),Downsample), ...
                                       movmean(CellData{m,1}(:),Downsample),NewTime(:),Method);
           else
               CellData{m,1} = NewTime;
               if strcmp(FieldNames{m,1},'TimeStamp') ~= true
                   CellData{m,1} = CellData{m,1}.*nan;
               end
           end
       end
       CellDataNew{m,1} = CellData{m,1};
   end
end
end




