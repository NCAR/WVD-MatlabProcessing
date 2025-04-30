% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created July 28th, 2020

function [Data,LidarData] = UnpackRawData(RawData)
%
% Inputs: RawData:    A data structure that contains all of the loaded data
%
% Outputs: Data:      A data structure containing all loaded data unpacked
%                     into a usable form
%          LidarData: A structure breaking out the photon count data
%                     produced by the lidar (not life support data)
%
%% Defining elements to unpack
ToUnpack  = {'Container';'Etalon';'Laser';'MCS';'Current';'Thermocouple'};
MCSField  = 'MCS';
SubField  = 'Data';
%% Unpacking power data first
[IsField,Power] = RecursivelyCheckIsField(RawData,'Power');
if IsField
    % Parsing out channel types for MCX...not necessary for the MCS
    if size(Power.LaserPower,2) ~= size(Power.Type,2)
        % Figuring out which block this power data belongs to
        ChannelNum = size(Power.LaserPower,2);
        Start      = ChannelNum.*Power.MCSInstance + 1;
        % Pre-allocating data array
        Type2 = nan.*ones(size(Power.LaserPower,1),ChannelNum);
        % Looping over and removing pulse blocks
        for m=1:1:size(Power.LaserPower,1)
            Type2(m,:) = Power.Type(m,Start(m):(Start(m)+ChannelNum-1));
        end
        % Overwritting type with the actual value
        Power.Type = Type2;
    end
    % Pulling out the data label from the array
    Key = Power.Type; Key(isnan(Key)) = []; PTypes = unique(Key);
    PTypes(PTypes==PowerChannels({'Unrecognized'})) = [];
    % Checking if the key is default data (nan)
    if isempty(Key)
        Key = ones(size(Power.Type)).*PowerChannels({'Unrecognized'});
    end
    % Converting power data to cell array
    Power  = rmfield(Power,'Type'); Power  = rmfield(Power,'MCSInstance');
    FN     = fieldnames(Power);
    Power  = struct2cell(Power);
    % Looping over all data and parsing
    for m=1:1:size(PTypes,1)     % Looping over power types
        % Determining if there is data from this data type in this row
        DataExistsInRow = sum(Key == PTypes(m),2) >= 1;
        for n=1:1:size(Power,1)  % Looping over fieldnames
            % Temporary variables to remove rows without relevant data
            TempVarKeyData   = Key(DataExistsInRow,:);
            TempVarPowerData = Power{n,1}(DataExistsInRow,:);
            % Saving data
            if all(size(TempVarKeyData) ==  size(TempVarPowerData))
                Temp{m,1}{n,1} = TempVarPowerData(TempVarKeyData == PTypes(m));
            else
                Temp{m,1}{n,1} = TempVarPowerData; %#ok<*AGROW>
            end
        end
    end
    % Converting to a full structure and cleaning up unneeded variables
    RawData.Power = cell2struct(cellfun(@(s) cell2struct(s,FN),Temp,'Uni',false),PowerChannels(PTypes));
    clear EqualKey FN IsField Key m n PTypes Temp
end
clear IsField Power
%% Converting structure to cell to be able to loop over elements
FieldNames = fieldnames(RawData);
Data       = struct2cell(RawData);
%% Unpacking the container/etalon/laser/MCS/Current/Thermocouple data
for m=1:1:size(ToUnpack,1)
    Member = find(ismember(FieldNames,ToUnpack{m,1}));
    if isempty(Member) ~= 1
        Data{Member,1} = UnpackTimeSeriesData(Data{Member,1});
    end
end
%% Converting back to strucutre
Data = cell2struct(Data,FieldNames);
%% Checking for MCS data and handling accordingly 
[IsField,LidarData.Raw] = RecursivelyCheckIsField(Data,MCSField);
if IsField
    Fields = fieldnames(LidarData.Raw);
    for m=1:1:length(Fields)
        Data.(MCSField).(Fields{m}) = rmfield(Data.(MCSField).(Fields{m}),SubField);
    end
end
end

% Unpacking time series data containing multiple unique tags
function [NewData] = UnpackTimeSeriesData(Data)
%
% Inputs: Data:     A structure containing hardware data to be unpacked
%
% Outputs: NewData: A reformatted structure containing unpacked data
%
%% Finding unique hardware types
UniqueTypes = sortrows(unique(Data.Type));
AllTypes    = Data.Type;
%% Checking if any data exists that is not just default
if isa(AllTypes,'double')
    AllTypes(isnan(AllTypes)) = [];
    if isempty(AllTypes)
        AllTypes = cell(size(Data.Type));
        AllTypes(:) = {'Unrecognized'};
        UniqueTypes = {'Unrecognized'};
    end
end
%% Removing types from the strucutre and making it a cell array for looping
Data       = rmfield(Data,'Type');
FieldNames = fieldnames(Data);
Data       = struct2cell(Data);
%% Looping over hardware types and cell array elements to parse by hardware
NewData = cell(size(UniqueTypes,1),1);
% Checking if file elements need to be broadcast to multple elements...this
% happens when something writes multiple values per time stamp (like
% current monitoring) 
if size(AllTypes,2) ~= 1
    for m=1:1:size(Data,2)
        if size(Data{m,1},2) == 1
            Data{m,1} = repmat(Data{m,1},1,size(AllTypes,2));
        end
    end
end
% Parsing data
for m=1:1:size(UniqueTypes,1)
    % Finding the indices of each piece of hardware
    IsPresent = cellfun(@(s) strcmp(UniqueTypes{m,1},s), AllTypes);
    % Parsing original cell array into hardware specific sub-cell arrays
    if size(AllTypes,2) ~= 1
        NewData{m,1} = cellfun(@(s) s(IsPresent),Data,'Uni',false);
    else
        NewData{m,1} = cellfun(@(s) s(IsPresent,:),Data,'Uni',false);
    end
    % Converting the sub-cells back to structures
    NewData{m,1} = cell2struct(NewData{m,1},FieldNames);
end
%% Converting to a full structure
NewData = cell2struct(NewData,UniqueTypes);
end
