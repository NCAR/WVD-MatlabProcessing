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
Bad = {'Current'       ,'TimeStamp';
       'Etalon'        ,'TimeStamp';
       'HumiditySensor','TimeStamp';
       'Laser'         ,'TimeStamp';
       'Power'         ,'TimeStamp';
       'Thermocouple'  ,'TimeStamp';
       'UPS'           ,'TimeStamp';
       'WeatherStation','TimeStamp'};
% Don't currently see: Quantum Composer, Container, MCS
%%
for m=1:1:size(Bad,1)
    % Checking if the fields exist
    [IsField,OldTime] = RecursivelyCheckIsField(Data, Bad(m,:));
    % Field exists so check if it is monotonically increasing in time
    if IsField
       Data.(Bad{m,1}).(Bad{m,2}) = ForceMonotonicTimeStamps(OldTime);
       % Checking to see if the data needs to be shifted backwards
       for n=1:1:5  % Only moving back 5 days (more than 1 would be weird)
           if mean(Data.(Bad{m,1}).(Bad{m,2}),'omitnan') > 24 % Shifting back a full day to check
               Data.(Bad{m,1}).(Bad{m,2}) = Data.(Bad{m,1}).(Bad{m,2}) - 24;
           else
               break
           end
       end
       
    end
end
end

