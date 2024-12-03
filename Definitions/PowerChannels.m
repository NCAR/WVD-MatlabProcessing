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
S = struct('WVOffline',1,'WVOnline',2,'HSRL',3,'O2Offline',4,'O2Online',5,'Unrecognized',6);
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