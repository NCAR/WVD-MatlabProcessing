% Written by: Robert Stillwell 
% Written for: NCAR
% Modificication Info: Created February 13, 2018

function [ReturnDataStructure] = RecursivelyInterpolate2DStructure(OriginalDataStructure,NewTime,NewRange,Method)
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
[Cell,FieldNames,TimeStamps,Range] = RecursiveStruct2Cell(OriginalDataStructure);
%% Recursively performing an interpolation of the cell contents
CellDataNew = RecursiveInterpolateData(Cell,TimeStamps,NewTime,Range,NewRange,Method,FieldNames);
%% Convert the surface weather cell array back to a structure
ReturnDataStructure = RecursiveCell2Struct(CellDataNew, FieldNames);
end

% This function performs interpolation on all elements in a cell array
function [CellDataNew] = RecursiveInterpolateData(CellData,OldTime,NewTime,OldRange,NewRange,Method,FieldNames) 
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
       Temp = RecursiveInterpolateData(CellData{m,1},OldTime{m,1},NewTime,OldRange{m,1},NewRange,Method,FieldNames{m,2});
       CellDataNew{m,1} = Temp;                       %#ok<*AGROW>
   else
%        fprintf(['     ',FieldNames{m,1},'\n'])
       % The bottom of the tree so interpolate if it's not a single number
       if strcmp(FieldNames{m,1},'TimeStamp')
           CellData{m,1} = NewTime;
       elseif strcmp(FieldNames{m,1},'Range')
           CellData{m,1} = NewRange;
       else
           if all([size(OldTime(:),1),size(OldRange(:),1)] == size(CellData{m,1}))
               % Data needs to be transposed
               CellData{m,1} = CellData{m,1}';
           elseif all([size(OldRange(:),1),size(OldTime(:),1)] == size(CellData{m,1})) ~= 1
               % Data size is weird
               12
           end
           % Downsample data in time
           A = OldTime(:); B = NewTime(:);
           Downsample = ceil((B(2)-B(1))/(A(2)-A(1)));
           OldTime = movmean(A,Downsample,1) - Downsample*(A(2)-A(1));
           CellData{m,1} = movmean(CellData{m,1},Downsample,1);
           % Downsample data in range
           A = OldRange(:); B = NewRange(:);
           Downsample = floor((B(2)-B(1))/(A(2)-A(1)));
           OldRange = movmean(A,Downsample) - Downsample*(A(2)-A(1));
           CellData{m,1} = movmean(CellData{m,1},Downsample,2);
           % Interpolating in 2D
           [C,D] = meshgrid(NewTime,NewRange);
           CellData{m,1} = interp2(OldTime,OldRange,CellData{m,1},C,D);
       end
       CellDataNew{m,1} = CellData{m,1};
   end
end
end
