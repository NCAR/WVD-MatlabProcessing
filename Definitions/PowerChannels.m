% Written By: Robert Stillwell
% Written For: NCAR
% 
% Function to convert power channel map to and from integrers. This is done
% because the search for integers is much much faster than comparing
% strings.
function [Converted] = PowerChannels(Data)
%
% 
%
%% Defining the numeric codes
A = {'WVOffline'   ,'WVOnline'   ,'HSRL'   ,'O2Offline'   ,'O2Online',   ...
     'WVOfflineTWA','WVOnlineTWA','HSRLTWA','O2OfflineTWA','O2OnlineTWA'};
B = {'Unrecognized','Unassigned'};
C = arrayfun(@(m) sprintf('Unknown%02.0f',m),1:99,'UniformOutput',false);
% Converting to structure
TypeArray = [(1:length(A)),ones(1,length(B)).*length(A)+1,ones(1,length(C)).*length(A)+1];
S = cell2struct(num2cell(TypeArray),cat(2,A,B,C),2);

%% Converting
if isa(Data,'cell')  % Converting from string to number
    try % Convert all data strings to numbers 
        Converted = cellfun(@(s) S.(s),Data);
    catch
        CWLogging('****PowerChannels unrecognized****\n',Options,'Warning')
        % If there's a string that is unrecognized, remove it and try again
        FieldNames = fieldnames(S);
        Data(cellfun(@(s) sum(contains(FieldNames,s)),Data)==0) = [];
        Converted = cellfun(@(s) S.(s),Data);
    end
elseif isa(Data,'double') % Converting from number to string
    FieldNames = fieldnames(S);
    Data(Data>size(FieldNames,1)) = S.Unrecognized;
    Converted  = FieldNames(Data);
end
end