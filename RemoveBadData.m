% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August 20th, 2020

function [Data] = RemoveBadData(Data)
%
% Inputs: Data:  A data structure that contains all of the loaded data
%
% Outputs: Data: A data structure that has flagged bad data removed
%
%% Defining elements to remove (Top level, Sub Level, Test) 
Bad = {'Laser'         ,'WavelengthActual',-1;
       'HumiditySensor','RelativeHumidity',-100;
       'Thermocouple'  ,'Temperature'     ,-1000; % Broadcast error
       'WeatherStation','RelativeHumidity',-1e-9;  
       'Container'     ,'TimeStamp'       , 1e-9};
%% Searching to see if the elements of the Bad data array exist
for m=1:1:size(Bad,1)
    % Checking if the fields exist
    [IsField,Data2Check] = RecursivelyCheckIsField(Data,Bad(m,1:2));
    % Field exists so check if there's bad data...remove as needed
    if IsField
        Data.(Bad{m,1}) = RemoveData(Data.(Bad{m,1}),Data2Check>Bad{m,3});
    end
end    
end

% This function removes identified bad data from the data structure
function [Struct] = RemoveData(Struct,GoodData)
%
% Inputs: Struct:   A structure containing raw data read from the MPD files
%         GoodData: An array of boolean objects indicating good data in the 
%                   data structure to be preserved
%
% Outputs: Struct: A structure containing raw data that not marked as bad 
%
%% Checking if all data is bad (making sure not to return an empty struct)
if sum(~GoodData) == length(GoodData)
    [IsField,TimeStamp] = RecursivelyCheckIsField(Struct,'TimeStamp');
    % Removing bad data in the entire structure
    Cell = cellfun(@(X) X.*nan,struct2cell(Struct),'Uni',false);
    % Converting back to a named structure from a cell arrray
    Struct = cell2struct(Cell,fieldnames(Struct));
    % Putting back the original timestamps
    if IsField
        Struct.TimeStamp = TimeStamp;
    end
    return
end
%% Removing only the bad data
if sum(~GoodData) > 0
    % Checking if data has a broadcast time array
    [~,C] = cellfun(@size,struct2cell(Struct),'UniformOutput',true);
    % If all elements aren't the same size, assume error can be broadcast
    if ~all(C(1)==C)
        GoodData = GoodData(:,1);
    end
    % Removing bad data in the entire structure
    Cell = cellfun(@(X) X(GoodData,:),struct2cell(Struct),'Uni',false);
    % Converting back to a named structure from a cell arrray
    Struct = cell2struct(Cell,fieldnames(Struct));
end
end
