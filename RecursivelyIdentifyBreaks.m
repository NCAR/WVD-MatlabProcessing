% Written by: Robert Stillwell 
% Written for: NCAR
% Modificication Info: Created August 1, 2018

function [ReturnDataStructure] = RecursivelyIdentifyBreaks(OriginalDataStructure,BreakSize)
%
% Inputs: OriginalDataStructure: Data structure to search for data breaks
%         BreakSize:             Medians allowed before marking data break
%
% Outputs: ReturnDataStructure:  A structure containing all data as before
%                                but with data breaks padded with NaNs
%
%% Converting the surface weather structure into a cell array
[Cell,FieldNames,TimeStamps] = RecursiveStruct2Cell(OriginalDataStructure);
%% Recursively performing an interpolation of the cell contents
CellDataNew = RecursiveFillBreaks(Cell,TimeStamps,BreakSize);
%% Convert the surface weather cell array back to a structure
ReturnDataStructure = RecursiveCell2Struct(CellDataNew, FieldNames);
end

% This function fills data breaks with a single NaN so it won't interpolate
% over breaks
function [CellDataNew] = RecursiveFillBreaks(CellData, OldTime, BreakSize) 
%                        
% Inputs: CellData:      Cell array to recursively search for data breaks
%         OldTime:       Cell array of time stamps from each cell level
%         BreakSize:     Medians allowed before marking data break
%                        
% Outputs: CellDataNew:  Cell array with all data as before but with data 
%                        breaks padded with NaNs
%                        
%% Recursively checking data contained within the cell array for breaks
for m=1:1:size(CellData,1)
   if iscell(CellData{m,1}) 
       % Need to dive down further into the cell array
       Temp = RecursiveFillBreaks(CellData{m,1},OldTime{m,1},BreakSize);
       CellDataNew{m,1} = Temp; %#ok<*AGROW>
   else
       % At the bottom of the cell tree so check data for breaks
       if size(OldTime,1) == size(CellData{m,1},1)
           % Identifying breaks larger than the allowable 
           Differences = diff(OldTime);
           MedDiff     = median(Differences);
           Breaks      = find(Differences > MedDiff.*BreakSize);
           if isempty(Breaks) == 0
               % Breaks found so fill them
               IsTimeArray = all(OldTime == CellData{m,1}(:,1));
               for n=1:1:length(Breaks)
                   Element2FIll = nan(1,size(CellData{m,1},2));
                   if IsTimeArray % Array is a timestamp so can't fill with nan
                       Element2FIll = ones(size(Element2FIll)).* ...
                                      CellData{m,1}(Breaks(n)+n-1)+MedDiff;
                   end
                   CellData{m,1} = [CellData{m,1}(1:(Breaks(n)+n-1),:);Element2FIll;
                                    CellData{m,1}((Breaks(n)+n):end,:)];
               end
           end
       end
       CellDataNew{m,1} = CellData{m,1};
   end
end
end
