% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created March 19, 2020

function [Data,Found] = ReadMPDData(DataPath,Options)
%
% Input: DataPath: String containing the full path of the data
%        Options:  A structure containing all user defined options desired
%                  for processing MPD data
%                  
% Output: Data:    Structure containing all of the requested data if found
%
%% Defining the structure of the data information cell array
CType.File     = 1; CType.DataName = 2; CType.FileVar  = 3;
CType.CodeVar  = 4; CType.VarType  = 5; Found = false;
%% Determining the file structure
ToLoad = DefineFileStructure(Options.DataNames,CType);
%% Loading data
for m=1:1:size(ToLoad,1)                              % Looping over filetypes
    s = dir(fullfile(DataPath,ToLoad{m,CType.File})); % Finding availible files
    % Loading all files
    if isempty(s)
        % No files to load
        TimeBounds = linspace(0,24,100)';
        % Looping over varaibles and filling with default blank data
        for p=1:1:size(ToLoad{m,CType.FileVar})
            CWLogging(['Loading: ',ToLoad{m,CType.File},', default\n'],Options,'Sub')
            % Pre-allocating data array space
            if strcmp(ToLoad{m,CType.CodeVar}{p},'TimeStamp')
                Data{m,1}{p,1} = TimeBounds;
            else
                Data{m,1}{p,1} = TimeBounds.*nan; %#ok<*AGROW>
            end
        end
    else
        Found = true;
        % Looping over varaibles
        for p=1:1:size(ToLoad{m,CType.FileVar})
            % Writing logging information to the command window
            CWLogging(['Loading: ',ToLoad{m,CType.File},', ',ToLoad{m,CType.FileVar}{p},...
                           ' as ',ToLoad{m,CType.CodeVar}{p},'\n'],Options,'Sub')
            % Looping over availible files and loading data
            for n=1:1:size(s,1)
                Filename     = fullfile(DataPath,s(n,1).name);
                FileVar      = ToLoad{m,CType.FileVar}{p};
                FileType     = ToLoad{m,CType.DataName};
                VariableType = ToLoad{m,CType.VarType}{p};
                % Reading data file
                A = ReadVariable(Filename,FileVar,VariableType);
                % Check if there are special instructions with this load
                TempData{n,1} = SpecialInsructions(A, Filename, FileVar,FileType);
            end
            % Checking data is all the same size (or padding) & combining
            [~,DataBins] = cellfun(@size,TempData);
            MaxDim = max(DataBins);
            if any(DataBins ~= MaxDim)
                CWLogging('      Warning: Resizing data\n',Options,'Sub')
                % Padding the end of the array with nan vallues
                TempData = cellfun(@(X) [X(:,:),nan(size(X,1),MaxDim-size(X,2))],TempData,'Uni',false);
            end
            Data{m,1}{p,1} = vertcat(TempData{:});
            clear TempData
        end  
    end
end
%% Converting Raw data cell array to a strucutre
for m=1:1:size(Data,1)
    Data{m,1} = cell2struct(Data{m,1},ToLoad{m,CType.CodeVar});
end
Data = cell2struct(Data,Options.DataNames);
end

% This function is used to simply define the structure of all files needed
% for MPD processing
function [ToLoad] = DefineFileStructure(Types,CType)
%
% Inputs: Types:   A cell array containing the file name types the user
%                  desires to load
%         CTypes:  Definitions of the column types for the file structure 
%                  cell array
%
% Outputs: ToLoad: A cell array containing information on all possible data
%                  types and how the files that contain them are formatted
%
%% Looping over the file types and determining their structure 
for m=1:1:length(Types)
    FileType = Types{m};
    % Defining the structure of the netcdf file
    switch FileType
        case 'Container'
            DataType = 'Container*.nc';
            Vars     = {'time';'FunctionType';'CQueueEl';'RQueueEl';'Status'};
            CodeName = {'TimeStamp';'Type';'-';'-';'-'};
            Type     = {'-';'String';'-';'-';'-'};
        case 'Current'
            DataType = 'Current*.nc';
            Vars     = {'time';'Current';'CurrentMonitorLocations'};
            CodeName = {'TimeStamp';'-';'Type'};
            Type     = {'-';'-';'String'};
        case 'Etalon'
            DataType = 'Etalon*.nc';
            Vars     = {'time';'Temperature';'TempDiff';'IsLocked';'EtalonNum'};
            CodeName = {'TimeStamp';'TemperatureActual';'-';'-';'Type'};
            Type     = {'-';'-';'-';'-';'String'};
        case 'HumiditySensor'
            DataType = 'Humidity*.nc';
            Vars     = {'time';'InternalTemperature';'ExternalTemperature';
                        'RelativeHumidity'};
            CodeName = {'TimeStamp';'IntTemp';'ExtTemp';'-'};
            Type     = {'-';'-';'-';'-'};
        case 'Laser'
            DataType = 'LL*.nc';
            Vars     = {'time';'Wavelength';'WaveDiff';'IsLocked';'TempDesired';
                        'TempMeas';'Current';'SeedPower';'LaserName'};
            CodeName = {'TimeStamp';'WavelengthActual';'-';'Locked';'TemperatureDesired';
                        'TemperatureActual';'-';'-';'Type'};
            Type     = {'-';'-';'-';'-';'-';'-';'-';'-';'String'};
        case 'MCS'
            DataType = 'MCS*.nc';
            Vars     = {'time';'ProfilesPerHist';'Channel';'nsPerBin';'NBins';
                        'Data';'RTime'};
            CodeName = {'TimeStamp';'ProfilesPerHistogram';'Type';'RangeResolution';
                        '-';'-';'-'};
            Type     = {'-';'-';'-';'-';'-';'-';'-'};
        case 'Power'
            DataType = 'Power*.nc';
            Vars     = {'time';'RTime';'Power';'AccumEx';'Demux';'ChannelAssignment'};
            CodeName = {'TimeStamp';'-';'LaserPower';'-';'-';'Type'};
            Type     = {'-';'-';'-';'-';'-';'String'};
        case 'QuantumComposer'
            DataType = 'Clock*.nc';
            Vars     = {'time';'PulseDelay';'GateDelay';'DutyCycle';'SwitchRate';'PulseDuration';'PRF';'RiseTime';
                        'TSOA';'Online';'Offline';'Gate'};
            CodeName = {'TimeStamp';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-'};
            Type     = {'-';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-'};
        case 'Thermocouple'
            DataType = 'HKeep*.nc';
            Vars     = {'time';'Temperature';'ThermocoupleLocations'};
            CodeName = {'TimeStamp';'-';'Type'};
            Type     = {'-';'-';'String'};
        case 'UPS'
            DataType = 'UPS*.nc';
            Vars     = {'time';'BatteryNominal';'BatteryReplace';'BatteryInUse';
                        'BatteryLow';'BatteryCapacity';'BatteryTimeLeft';
                        'UPSTemperature';'HoursOnBattery'};
            CodeName = {'TimeStamp';'-';'-';'-';'-';'-';'-';'Temperature';'-'};
            Type     = {'-';'-';'-';'-';'-';'-';'-';'-';'-'};
        case 'WeatherStation'
            DataType = 'WS*.nc';
            Vars     = {'time';'Temperature';'RelHum';'Pressure';'AbsHum'};
            CodeName = {'TimeStamp';'-';'RelativeHumidity';'-';'AbsoluteHumidity'};
            Type     = {'-';'-';'-';'-';'-'};
        otherwise
            DataType = '';
            Vars     = {};
            CodeName = {};
            Type     = {};
    end
    % Populating cells with "-" variables with their partner
    CodeName(contains(CodeName,'-')==1) = Vars(contains(CodeName,'-')==1);
    Type(contains(Type,'-')==1) = {'Double'};
    % Storing data 
    ToLoad{m,CType.File}     = DataType; %#ok<*AGROW>
    ToLoad{m,CType.FileVar}  = Vars;
    ToLoad{m,CType.CodeVar}  = CodeName;
    ToLoad{m,CType.VarType}  = Type;
    ToLoad{m,CType.DataName} = FileType;
end
end

% This function is used to define when a simple data load is not
% acceptable. This is used primarily when a second data variable is
% requried to interpet the first. 
function [Data] = SpecialInsructions(Data,Filename,FileVar,FileType)
%
% Inputs: Data:     Loaded data from a netcdf file
%         Filename: String containing the desired file name
%         FileVar:  NetCDF variable name to load
%         FileType: Type of file being loaded
%
% Outputs: Data:    Modified loaded data per the belwo instructions
%
%%
if strcmp(FileType,'Current') && strcmp(FileVar,'Current')
    % Current is written in the file is half the actual value
    Data = Data.*2;
end
%% MCS Channel (also load the channel map and convert raw data with it)
if strcmp(FileType,'MCS') && strcmp(FileVar,'Channel')
    % Loading the channel map
    Key = ReadVariable(Filename,'ChannelAssignment','String');
    % If Key is read without error, use it to convert the channel number
    if iscell(Key)
        % Pre-allocating cell array
        NewData = cell(length(Data),1);
        % Filling the cell array with key data
        for m=1:1:length(NewData)
            NewData{m,1} =  Key{Data(m,1)+1};
        end 
        % Returning key data instead of channel numbers
        Data = NewData;
    end
end
%% Power Channel Map (also load channel map to see what should be kept)
if strcmp(FileType,'Power')
    % Loading the channel map
    Key = ReadVariable(Filename,'ChannelAssignment','String');
    % If Key is read without error, use it to convert the channel number
    if iscell(Key) && length(Key) == size(Data,2)     % Power data  
        Data(:,contains(Key,'Unassigned')) = [];
    elseif iscell(Key) && length(Key) == size(Data,1) % Labeling data
        Data(contains(Key,'Unassigned')) = [];
        Data = repmat(PowerChannels(Data)',ncinfo(Filename,'time').Size,1);
    end
end
%% Reshaping the location variable 
if strcmp(FileType,'Thermocouple') && strcmp(FileVar,'ThermocoupleLocations')
    % Channel map is already loaded so try to load the number of
    % measurements to repmat the map over
    Key = ReadVariable(Filename,'Temperature','Double');
    if iscell(Data) ~= 1 || size(Key,2) ~= length(Data)
        for m=1:1:size(Key,2); DKey{m,1} = sprintf('Unknown%02.0f',m); end
        Data = DKey; % Setting data to a default value
    end
    Data = repmat(Data',size(Key,1),1);
end
if strcmp(FileType,'Current') && strcmp(FileVar,'CurrentMonitorLocations')
    % Channel map is already loaded so try to load the number of
    % measurements to repmat the map over
    Key = ReadVariable(Filename,'Current','Double');
    if iscell(Data) ~= 1 || size(Key,2) ~= length(Data)
        for m=1:1:size(Key,2); DKey{m,1} = sprintf('Unknown%02.0f',m); end
        Data = DKey; % Setting data to a default value
    end
    Data = repmat(Data',size(Key,1),1);
end
end

% This function tries to read a variable of interest. If the type is Double
% it will typecast the variable, otherwise the type is maintained from the 
% netcdf file. Finally, if the variable type is written as a string, h5read
% is required because ncread does not recognize the string type. 
function [A] = ReadVariable(Filename,FileVar,VariableType)
%
% Inputs: Filename:     String containing the desired file name
%         FileVar:      NetCDF variable name to load
%         VariableType: Data type of variable to read
%
% Outputs: A:           Loaded data
%
%% Checking data inputs 
if nargin == 2
    VariableType = 'Double';
end
%% Loading data
A = [];
try
    if strcmp(VariableType,'Double')
        A = double(ncread(Filename,FileVar));
    elseif strcmp(VariableType,'Native')
        A = ncread(Filename,FileVar);
    elseif strcmp(VariableType,'String')
        A = h5read(Filename,['/',FileVar]);
        % Matlab 2017 returns A as a cell array but matlab 2020 returns it 
        % as a string array. Downstream needs a cell.
        if isa(A,'string') 
            A = cellstr(A);  
        end
    end
catch
    fprintf('     Failed to load variable\n')
    A = nan;
end
end
