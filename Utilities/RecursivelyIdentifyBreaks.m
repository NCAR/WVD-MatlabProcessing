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
function [DataNew] = RecursiveFillBreaks(Data, OldTime, BreakSize)
%                        
% Inputs: Data:          Cell array to recursively search for data breaks
%         OldTime:       Cell array of time stamps from each cell level
%         BreakSize:     Medians allowed before marking data break
%                        
% Outputs: DataNew:      Cell array with all data as before but with data
%                        breaks padded with NaNs
%                        
%% Recursively checking data contained within the cell array for breaks
for m=1:1:size(Data,1)
   if iscell(Data{m,1})
       % Need to dive down further into the cell array
       DataNew{m,1} = RecursiveFillBreaks(Data{m,1},OldTime{m,1},BreakSize);
   else
       % At the bottom of the cell tree so check data for breaks
       if size(OldTime,1) == size(Data{m,1},1)
           % Identifying breaks larger than the allowable 
           Differences = diff(OldTime);
           MedDiff     = median(Differences);
           % Checking the break look time is not less than 2 minutes
           if (MedDiff.*BreakSize)*60 > 2
               BreakAllowed = MedDiff.*BreakSize;
           else
               BreakAllowed = 2/60; % 2 minute in units of hours
           end
           B           = find(Differences > BreakAllowed);
           if isempty(B) == 0
               % Breaks found so fill them
               IsTimeArray = all(OldTime == Data{m,1}(:,1));
               for n=1:1:length(B)
                   ToFill = nan(1,size(Data{m,1},2));
                   if IsTimeArray % Array is timestamp so can't fill with nan
                       ToFill = ones(size(ToFill)).*Data{m,1}(B(n)+n-1)+MedDiff;
                   end
                   Data{m,1} = [Data{m,1}(1:(B(n)+n-1),:);ToFill;Data{m,1}((B(n)+n):end,:)];
               end
           end
       end
       DataNew{m,1} = Data{m,1}; %#ok<*AGROW>
   end
end
end
