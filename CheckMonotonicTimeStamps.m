% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August 19th, 2020

function [Data] = CheckMonotonicTimeStamps(Data)
%
% Inputs: Data:  A data structure that contains all of the loaded data
%
% Outputs: Data: A data structure that has non-monitonically increasing
%                time stamps addressed.
%
%% Fields to check if they are monotonically increasing 
% Bad = {'MCS'       ,'TimeStamp';
%        'Power'     ,'TimeStamp';
%        'Container' ,'TimeStamp'};
Bad = {'Power'     ,'TimeStamp'};
%%
for m=1:1:size(Bad,1)
    % Checking if the fields exist
    [IsField,OldTime] = RecursivelyCheckIsField(Data, Bad(m,:));
    % Field exists so check if it is monotonically increasing in time
    if IsField
       Data.(Bad{m,1}).(Bad{m,2}) = ForceMonotonicTimeStamps(OldTime);
    end
end
end

