% Written by: Robert Stillwell
% Modified for: National Center For Atmospheric Research
% Modification info: Created: November 8, 2017


function [Counts,PulseInfoNew] = ReadRawNetCDFData(DataTypes,HardwareMap,Paths)
%
%
%
%
%
%
%% Printing status to the command window
fprintf('Loading Data\n')

%%
cd(Paths.RawNetCDFData)

%% Possible hardware types
EtalonTypes = {'WVEtalon','HSRLEtalon','O2Etalon'};
LaserTypes  = {'WVOnline','WVOffline','HSRL','O2Online','O2Offline'};

%% Pre-allocating data
Etalon.TemperatureActual  = []; Etalon.TemperatureDesired = [];
Etalon.TimeStamp          = []; Etalon.Type               = [];

Thermocouple.TimeStamp    = []; Thermocouple.Temperature  = [];

Laser.Current             = []; Laser.Locked              = [];
Laser.TemperatureActual   = []; Laser.TemperatureDesired  = [];
Laser.TimeStamp           = []; Laser.Type                = [];
Laser.WavelengthActual    = []; Laser.WavelengthDesired   = [];

MCS.Channel               = []; MCS.Data                  = [];
MCS.ProfilesPerHistogram  = []; MCS.RangeResolution       = []; 
MCS.TimeStamp             = [];
         
Power.LaserPower          = []; Power.TimeStamp           = [];

UPS.BatteryCapacity       = []; UPS.BatteryHours          = [];
UPS.BatteryInUse          = []; UPS.BatteryNominal        = [];
UPS.Temperature           = []; UPS.TimeStamp             = [];

WStation.AbsoluteHumidity = []; WStation.Pressure         = [];
WStation.RelativeHumidity = []; WStation.Temperature      = [];
WStation.TimeStamp        = [];

%% Loading data
for m=1:1:size(DataTypes,1)
    % Finding all of the relevant files
   s = dir(DataTypes{m,1}); 
   for n=1:1:size(s,1)
       % Determining the current file name
       Filename = s(n,1).name;
       % Loading the data
       switch DataTypes{m,1}
           case 'Etalonsample*.nc'
               A = h5read(Filename,'/EtalonNum');
               Type = -1.*ones(size(A));
               for p=1:1:size(A,1)   % Looping over all data entries
                   for q=1:1:size(EtalonTypes,2) % Looping over all types
                       if isequal(A{p,1},EtalonTypes{q})
                           % If the type is recognized, save it
                           Type(p,1) = q - 1;
                           break
                       end
                   end
               end
               Etalon.TimeStamp          = [Etalon.TimeStamp         ;double(ncread(Filename,'time'))];
               Etalon.Type               = [Etalon.Type              ;Type];
               Etalon.TemperatureActual  = [Etalon.TemperatureActual ;ncread(Filename,'Temperature')];
               Etalon.TemperatureDesired = [Etalon.TemperatureDesired;ncread(Filename,'TempDiff')];
               clear Type
           case 'HKeepsample*.nc'
               Thermocouple.TimeStamp   = [Thermocouple.TimeStamp;   double(ncread(Filename,'time'))];
               Thermocouple.Temperature = [Thermocouple.Temperature; double(ncread(Filename,'Temperature'))];
           case 'LLsample*.nc'
               A = h5read(Filename,'/LaserName');
               Type = -1.*ones(size(A));
               for p=1:1:size(A,1)   % Looping over all data entries
                   for q=1:1:size(LaserTypes,2) % Looping over all types
                       if isequal(A{p,1},LaserTypes{q})
                           % If the type is recognized, save it
                           Type(p,1) = q - 1;
                           break
                       end
                   end
               end
               Laser.TimeStamp           = [Laser.TimeStamp          ;double(ncread(Filename,'time'))];
               Laser.Type                = [Laser.Type               ;Type];
               Laser.WavelengthActual    = [Laser.WavelengthActual   ;double(ncread(Filename,'Wavelength'))];
               Laser.WavelengthDesired   = [Laser.WavelengthDesired  ;double(ncread(Filename,'WaveDiff'))];
               Laser.Locked              = [Laser.Locked             ;double(ncread(Filename,'IsLocked'))];
               Laser.TemperatureActual   = [Laser.TemperatureActual  ;ncread(Filename,'TempMeas')];
               Laser.TemperatureDesired  = [Laser.TemperatureDesired ;ncread(Filename,'TempDesired')];
               Laser.Current             = [Laser.Current            ;ncread(Filename,'Current')];
               clear Type
           case 'MCSsample*.nc'
               A = double(ncread(Filename,'time'));
               if str2double(Filename(10:11)) == 23
                   % Checking to make sure all data is reported today
                   A(A<23) = A(A<23)+24;
               end
               MCS.TimeStamp             = [MCS.TimeStamp            ;A];
               MCS.ProfilesPerHistogram  = [MCS.ProfilesPerHistogram ;ncread(Filename,'ProfilesPerHist')];
               MCS.Channel               = [MCS.Channel              ;ncread(Filename,'Channel')];
               try
                   MCS.RangeResolution   = [MCS.RangeResolution      ;ncread(Filename,'CntsPerBin')];
               catch
                   MCS.RangeResolution   = [MCS.RangeResolution      ;ncread(Filename,'nsPerBin')];
               end
               MCS.Data                  = [MCS.Data                 ;single(ncread(Filename,'Data'))];  
               clear A
           case 'Powsample*.nc'
               Power.LaserPower          = [Power.LaserPower         ;ncread(Filename,'Power')];
               A = double(ncread(Filename,'time'));
               if str2double(Filename(10:11)) == 23
                   % Checking to make sure all data is reported today
                   A(A<23) = A(A<23)+24;
               end
               Power.TimeStamp           = [Power.TimeStamp          ;A];
               clear A
           case 'UPSsample*.nc'
               UPS.BatteryCapacity       = [UPS.BatteryCapacity ; ncread(Filename,'BatteryCapacity')];
               UPS.BatteryHours          = [UPS.BatteryHours    ; ncread(Filename,'BatteryTimeLeft')];
               UPS.BatteryInUse          = [UPS.BatteryInUse    ; ncread(Filename,'BatteryInUse')];
               UPS.BatteryNominal        = [UPS.BatteryNominal  ; double(ncread(Filename,'BatteryNominal'))];
               UPS.Temperature           = [UPS.Temperature     ; ncread(Filename,'UPSTemperature')];
               UPS.TimeStamp             = [UPS.TimeStamp       ; ncread(Filename,'time')];
           case 'WSsample*.nc'
               A = double(ncread(Filename,'time'));
               if str2double(Filename(9:10)) == 00
                   A(A>23) = A(A>23)-24;
               end
               WStation.TimeStamp        = [WStation.TimeStamp       ;A];
               WStation.Temperature      = [WStation.Temperature     ;ncread(Filename,'Temperature')];
               WStation.RelativeHumidity = [WStation.RelativeHumidity;ncread(Filename,'RelHum')];
               WStation.Pressure         = [WStation.Pressure        ;ncread(Filename,'Pressure')];
               WStation.AbsoluteHumidity = [WStation.AbsoluteHumidity;ncread(Filename,'AbsHum')];
       end
   end
end
clear m n s A Filename
% Converting differences to actual desired numbers
Etalon.TemperatureDesired = Etalon.TemperatureActual - Etalon.TemperatureDesired;
Laser.WavelengthDesired   = Laser.WavelengthActual - Laser.WavelengthDesired;

%% Changing back to the coding directory
cd(Paths.Code)

%% Getting rid of identical time stamps in power monitoring data 
A = find(diff(Power.TimeStamp) == 0);
Counter = 0;
while ~isempty(A)
    Counter = Counter + 1;
    Power.TimeStamp(diff([0;Power.TimeStamp(1:end-1)]) == 0) = Power.TimeStamp(diff([0;Power.TimeStamp(1:end-1)]) == 0) + 1e-9;
    A = find(diff(Power.TimeStamp) == 0);
    if Counter == 300
        fprintf('More than 300 power time stamps were identical.\n')
        break
    end
end
clear A Counter

%% Getting rid of incomplete data
% Checking for MCS data
MCSFirst = find(MCS.Channel == min(HardwareMap.PhotonCounting),1,'first');
MCSLast  = find(MCS.Channel == max(HardwareMap.PhotonCounting),1,'last');
% Removing incomplete scan at the start of the day
if MCSFirst ~= 1
    MCS.Channel(1:(MCSFirst-1)) = nan;
end
% Removing incomplete scan at the end of the day
if MCSLast ~= size(MCS.Channel,1)
    MCS.Channel((MCSLast+1):end) = nan;
end
clear MCSFirst MCSLast
% Removing incomplete scan in the middle of the day
MCS = RemoveIncompleteMCSScans(MCS,1);

%% Downsampling power data to roughly MCS timegrid
% Determining how many data points to average
PowerMeasurementsPerData = ceil(size(Power.TimeStamp,1)./ ...
                                size(MCS.TimeStamp(MCS.Channel == MCS.Channel(1)),1));
% Downsampling the data by averaging then taking the center values
if size(Power.TimeStamp,1) > 0
%Power  = RecursivelyDownSample(Power,PowerMeasurementsPerData,size(Power.TimeStamp,1));
else
TimeBounds = linspace(0,24,100)';
Power.LaserPower = nan.*TimeBounds;
Power.TimeStamp  = TimeBounds;
end     

%% Removing bad laser scans
BadData = find(Laser.WavelengthActual <= -1);
if isempty(BadData) == 0
    % Converting the names structure to a cell array
    CellArray  = struct2cell(Laser);
    FieldNames = fieldnames(Laser);
    % Removing bad data
    for m=1:1:size(CellArray)
        CellArray{m,1}(BadData,:) = [];
    end
    % Converting back to a named structure
    Laser = cell2struct(CellArray,FieldNames);
    clear CellArray FieldNames m
end
clear BadData 

%% Filling in bad data
TimeBounds = linspace(0,24,100)';
if isempty(WStation.TimeStamp)
   WStation.TimeStamp        = TimeBounds;
   WStation.AbsoluteHumidity = TimeBounds.*nan;
   WStation.Pressure         = TimeBounds.*nan;
   WStation.RelativeHumidity = TimeBounds.*nan;
   WStation.Temperature      = TimeBounds.*nan;
end
if isempty(Thermocouple.TimeStamp)
   Thermocouple.TimeStamp    = TimeBounds;
   Thermocouple.Temperature  = TimeBounds.*nan;
end
if isempty(UPS.TimeStamp)
   UPS.TimeStamp             = TimeBounds;
   UPS.BatteryCapacity       = TimeBounds.*nan;
   UPS.BatteryHours          = TimeBounds.*nan;
   UPS.BatteryInUse          = TimeBounds.*nan;
   UPS.BatteryNominal        = TimeBounds.*nan;
   UPS.Temperature           = TimeBounds.*nan;
end

%% Just making sure to get rid of laser power off data or error returns 
% Especially DIAL01 seems to return current = 0 often...just remove data
% and move on
Laser.Current(Laser.Current == 0) = nan;

%% Marking time gaps
Laser        = PaddingDataStructureTimeSeries(Laser,5,0);
Thermocouple = PaddingDataStructureTimeSeries(Thermocouple,5,0);
Etalon       = PaddingDataStructureTimeSeries(Etalon,5,1);
WStation     = PaddingDataStructureTimeSeries(WStation,5,0);
UPS          = PaddingDataStructureTimeSeries(UPS,5,0);
Power        = PaddingDataStructureTimeSeries(Power,5,0);
MCS          = PaddingDataStructureMCS(MCS,5,7043);

%% Parsing data
[Counts,PulseInfoNew] = RawNetCDFDataParse(Etalon,Laser,MCS,Power,Thermocouple,UPS,WStation,HardwareMap);
%% Pushing data to a regular grid
[~,PulseInfoNew]      = RawNetCDFData2RegularGrid(PulseInfoNew);
end

function [Counts,PulseInfo] = RawNetCDFDataParse(Etalon,Laser,MCS,Power,Thermocouple,UPS,WStation,HardwareMap)
%
%
%
%
%
%
%%
for m=1:1:size(HardwareMap.ChannelName,1)
    % Parsing out lidar data meta-data
    PulseInfo.Data.ProfilesPerHistogram{m,1}= MCS.ProfilesPerHistogram(MCS.Channel == HardwareMap.PhotonCounting(m));
    PulseInfo.Data.RangeResolution{m,1}     = MCS.RangeResolution(MCS.Channel == HardwareMap.PhotonCounting(m));
    % Parsing out the etalon
    PulseInfo.Etalon.TemperatureActual{m,1}  = Etalon.TemperatureActual(Etalon.Type == HardwareMap.Etalon(m) | isnan(Etalon.Type));
    PulseInfo.Etalon.TemperatureDesired{m,1} = Etalon.TemperatureDesired(Etalon.Type == HardwareMap.Etalon(m) | isnan(Etalon.Type));
    PulseInfo.Etalon.TimeStamp{m,1}          = Etalon.TimeStamp(Etalon.Type == HardwareMap.Etalon(m) | isnan(Etalon.Type));
    % Parsing out the thermocouple data
    PulseInfo.Housekeeping.Temperature       = Thermocouple.Temperature;
    % Parsing out the laser
    PulseInfo.Laser.Current{m,1}            = Laser.Current(Laser.Type == HardwareMap.Laser(m) | isnan(Laser.Type));
    PulseInfo.Laser.Locked{m,1}             = Laser.Locked(Laser.Type == HardwareMap.Laser(m) | isnan(Laser.Type));
    PulseInfo.Laser.Power{m,1}              = Power.LaserPower(:,HardwareMap.Power(m)+1);
    PulseInfo.Laser.TemperatureActual{m,1}  = Laser.TemperatureActual(Laser.Type == HardwareMap.Laser(m) | isnan(Laser.Type));
    PulseInfo.Laser.TemperatureDesired{m,1} = Laser.TemperatureDesired(Laser.Type == HardwareMap.Laser(m) | isnan(Laser.Type));
    PulseInfo.Laser.WavelengthActual{m,1}   = Laser.WavelengthActual(Laser.Type == HardwareMap.Laser(m) | isnan(Laser.Type));
    PulseInfo.Laser.WavelengthDesired{m,1}  = Laser.WavelengthDesired(Laser.Type == HardwareMap.Laser(m) | isnan(Laser.Type));
    % Parsing out the lidar data
    Counts.Raw{m,1}                         = MCS.Data(MCS.Channel == HardwareMap.PhotonCounting(m),:);
    % Time stamps
    PulseInfo.TimeStamp.Etalon{m,1}         = Etalon.TimeStamp(Etalon.Type == HardwareMap.Etalon(m) | isnan(Etalon.Type));
    PulseInfo.TimeStamp.Housekeeping        = Thermocouple.TimeStamp;
    PulseInfo.TimeStamp.LidarData{m,1}      = MCS.TimeStamp(MCS.Channel == HardwareMap.PhotonCounting(m));
    PulseInfo.TimeStamp.LaserLocking{m,1}   = Laser.TimeStamp(Laser.Type == HardwareMap.Laser(m) | isnan(Laser.Type));
    PulseInfo.TimeStamp.LaserPower{m,1}     = Power.TimeStamp;
    PulseInfo.TimeStamp.UPS                 = UPS.TimeStamp;
    PulseInfo.TimeStamp.WeatherStation      = WStation.TimeStamp;
    % Parsing out the UPS data
    PulseInfo.UPS.BatteryCapacity = UPS.BatteryCapacity;
    PulseInfo.UPS.BatteryHours    = UPS.BatteryHours;
    PulseInfo.UPS.BatteryInUse    = UPS.BatteryInUse;
    PulseInfo.UPS.BatteryNominal  = UPS.BatteryNominal;
    PulseInfo.UPS.Temperature     = UPS.Temperature;
    % Parsing out the weather station
    PulseInfo.WeatherStation.AbsoluteHumidity  = WStation.AbsoluteHumidity;
    PulseInfo.WeatherStation.Pressure          = WStation.Pressure;
    PulseInfo.WeatherStation.RelativeHumidity  = WStation.RelativeHumidity;
    PulseInfo.WeatherStation.Temperature       = WStation.Temperature;
end
end

function [PulseInfoOld,PulseInfo] = RawNetCDFData2RegularGrid(PulseInfo)
%
%
%
%
%%
PulseInfoOld = PulseInfo;
PulseInfo.TimeStamp.Merged = double(PulseInfo.TimeStamp.LidarData{1,1});
% Pushing weather station, housekeeping, and UPS data to MCS time grid
PulseInfo = RecursivelyInterpolateStructure(PulseInfo,PulseInfo.TimeStamp.WeatherStation,PulseInfo.TimeStamp.Merged,'linear','extrap');
PulseInfo = RecursivelyInterpolateStructure(PulseInfo,PulseInfo.TimeStamp.Housekeeping,PulseInfo.TimeStamp.Merged,'linear','extrap');
PulseInfo = RecursivelyInterpolateStructure(PulseInfo,PulseInfo.TimeStamp.UPS,PulseInfo.TimeStamp.Merged,'linear','extrap');

% Pushing power data to MCS time grid
PulseInfo = RecursivelyInterpolateStructure(PulseInfo,PulseInfo.TimeStamp.LaserPower{1,1},PulseInfo.TimeStamp.Merged,'linear','extrap');
% Pushing laser locking data to MCS time grid by looping over lasers...need
% to loop because native grids are all unique until this step
for m=1:1:size(PulseInfo.TimeStamp.LaserLocking)
    if size(PulseInfo.TimeStamp.LaserLocking{m,1},1) ~= size(PulseInfo.TimeStamp.Merged,1)
        PulseInfo = RecursivelyInterpolateStructure(PulseInfo,PulseInfo.TimeStamp.LaserLocking{m,1},PulseInfo.TimeStamp.Merged,'linear','extrap');
    end
end
% Pushing laser locking data to MCS time grid by looping over etalons
for m=1:1:size(PulseInfo.TimeStamp.Etalon)
    if size(PulseInfo.TimeStamp.Etalon{m,1},1) ~= size(PulseInfo.TimeStamp.Merged,1)
        PulseInfo = RecursivelyInterpolateStructure(PulseInfo,PulseInfo.TimeStamp.Etalon{m,1},PulseInfo.TimeStamp.Merged,'linear','extrap');
    end
end
end

% This is a general function used to identify missing data times when the
% data are not sent in serial. It inserts nan values into the time series
% which can not be interpolated thus marking all data gaps
function [DataPadded] = PaddingDataStructureTimeSeries(Data,MedianWidths2Flag,IgnoreTight)
%
%
%
%
%
%
%% Constants
TimeStepToAdd = 1e-9;
%% Defining function handles (distance in medians from the median)
FindOutliers  = @(DiffTimes) (DiffTimes - median(DiffTimes))./median(DiffTimes);
FindOutliers2 = @(DiffTimes,Spacing) (DiffTimes - Spacing)./Spacing;
%% Finding out if closely spaced readings are to be ignored
if IgnoreTight
    % Finding the median spacing of unique elements
    DiffTimes = diff([0;Data.TimeStamp;24]);
    Spacing = median(diff([0;Data.TimeStamp(Data.Type==min(unique(Data.Type)));24]));
    % With the second most common value, checking how far away other points
    % are from that in units of that
    A = find(FindOutliers2(DiffTimes,Spacing) >= MedianWidths2Flag);
else
    % Finding the points where the time stamps are more than ## median 
    % widths away from the median
    A = find(FindOutliers(diff([0;Data.TimeStamp;24])) >= MedianWidths2Flag);
end
% %% Finding the points where the time stamps are more than 5 median widths away from the median
% A = find(FindOutliers(diff([0;Data.TimeStamp;24])) >= MedianWidths2Flag);
%% Converting the names structure to a cell array
CellArray  = struct2cell(Data);
FieldNames = fieldnames(Data);
%% Looping over cell elements and filling them (if necessary)
if isempty(A)
    NewCell = CellArray;
else
    % Pre-allocating new data structure
    NewCell    = cell(size(CellArray));
    % Looping over cells to fill them with empty data
    for m=1:1:size(CellArray,1)
        % Initializing the new data array
        NewCell{m,1} = [];
        % Initializing the array location
        StartIndexOld = 1;
        EndOfDayBreak = 0;
        % Looping over all data gaps
        for n=1:1:size(A,1)
            EndIndexOld = A(n)-1;
            %%%%%%%% Checking if time stamps or nans should be added %%%%%%%%
            if strcmp(FieldNames{m,1},'TimeStamp')
                ToAdd = TimeStepToAdd;
            else
                ToAdd = nan;
            end
            %%%%%%%% Inserting data into time gaps %%%%%%%%
            if A(n) == 1
                %%%%%%%% There is a break at the start of the day %%%%%%%%
                % Determining if there are multiple breaks in the day,
                % otherwise, just make the end index the last measurement
                if size(A,1)>1
                    EndIndexOld = A(n+1)-1;
                else
                    EndIndexOld = size(CellArray{m,1},1);
                end
                % Padding the new cell array
                NewCell{m,1} = [NewCell{m,1};
                                CellArray{m,1}(StartIndexOld,:)-ToAdd;
                                CellArray{m,1}(StartIndexOld:EndIndexOld,:)];
                StartIndexOld = EndIndexOld+1;
            elseif A(n) == size(CellArray{m,1},1) + 1
                %%%%%%%% There is a break at the end of the day %%%%%%%%
                NewCell{m,1} = [NewCell{m,1};
                                CellArray{m,1}(StartIndexOld:EndIndexOld,:);
                                CellArray{m,1}(EndIndexOld,:)+ToAdd];
                EndOfDayBreak = 1;
            else
                %%%%%%%% There is a break in the middle of the day %%%%%%%%
                % Insert meaningless time stamps in the middle of the day
                NewCell{m,1} = [NewCell{m,1};
                                CellArray{m,1}(StartIndexOld:EndIndexOld,:);
                                CellArray{m,1}(EndIndexOld,:)+ToAdd;
                                CellArray{m,1}(EndIndexOld+1,:)-ToAdd];
                StartIndexOld = EndIndexOld+1;
            end
            %%%%%%%% FIlling end of day if needed %%%%%%%%
            if n == size(A,1) && EndOfDayBreak == 0
                NewCell{m,1} = [NewCell{m,1};
                                CellArray{m,1}(StartIndexOld:end,:)];
            end
        end
    end
end
%% Converting back to a named structure
DataPadded = cell2struct(NewCell,FieldNames);
end

% This function looks at the time series of MCS data and attempts to find
% any sets of data thclose aat are incomplete. This can happen when the system is
% turned on or off in the middle of the transmission of a set of MCS
% data-grams or when the set is interupted by the start or end of a day
function [MCS] = RemoveIncompleteMCSScans(MCS,Depth)
%
%
%
%
%
%
%% Removing data to simulate bad scans
% for m=1:1:6
%     Index2Flip(m) = floor(rand.*size(GoodRows,1)); %#ok<SAGROW>
% end
% GoodRows(Index2Flip) = [];
%% Identifying bad data
% Checking for the unique channel identification numbers
A = unique(MCS.Channel);
A(isnan(A)) = [];
% Pre-allocating data
BadData = [];
% Looping over the 
for m=1:1:size(A,1)
    % Find all differences different than the scan size
    ThisChannel = find(MCS.Channel == A(m));
    Incomplete  = find(diff(ThisChannel) < size(A,1));
    if isempty(Incomplete) == 0
        BadData     = [BadData;ThisChannel(Incomplete)]; %#ok<AGROW>
    end
end
BadData = [BadData;find(isnan(MCS.Channel))];
%% Printing the number of removed data points to the command window
if isempty(BadData) == 0
   fprintf('%0.0f measurements were removed as incomplete scans.\n',size(BadData,1)) 
end
%% Converting the names structure to a cell array
CellArray  = struct2cell(MCS);
FieldNames = fieldnames(MCS);
%% Removing bad data
for m=1:1:size(CellArray)
    CellArray{m,1}(BadData,:) = [];
end
%% Converting back to a named structure
MCS = cell2struct(CellArray,FieldNames);

%% Recursively checking for extra incomplete scans
if Depth < 10 && isempty(BadData) == 0
    % Removing more incomplete scans
    MCS = RemoveIncompleteMCSScans(MCS,Depth + 1);
elseif Depth > 10 && isempty(BadData) == 0
    % Notifying the user that the MCS scans are really really incomplete
    fprintf('More than 10 sets of incomplete scans were observed. Thats a problem.\n')
end
end

% This function specifically pads MCS data with bad values to remove any
% data breaks. This is different from the general function as it must
% handle 2d data and also must fill the channel type correctly otherwise
% the inserted data is not later understood
function [MCS] = PaddingDataStructureMCS(MCS,MedianWidths2Flag,LaserRepRate)
%
%
%
%
%
%% Defining function handles (distance in medians from the median)
FindOutliers = @(DiffTimes) (DiffTimes - median(DiffTimes))./median(DiffTimes);
%% Finding unique channel numbers
A = unique(MCS.Channel);
%% Finding the time stamps of the starts of the data transfers
ScanStarts = MCS.TimeStamp(1:size(A,1):end);
%% Finding the points where the time stamps are more than 5 median widths away from the median
B = find(FindOutliers(diff([0;ScanStarts;24])) >= MedianWidths2Flag);
DataBreakIndices = zeros(size(B,1),1);
for m=1:1:size(B,1)
    if B(m) <= size(ScanStarts,1) && B(m) ~= 1
        DataBreakIndices(m) = find(MCS.TimeStamp == ScanStarts(B(m)),1,'first');
    elseif B(m) == 1
        DataBreakIndices(m) = 1;
    else
        % Missing the end of the day
        DataBreakIndices(m)  = size(MCS.TimeStamp,1)+1;
        MCS.TimeStamp(end+1) = 24;
    end
end
%% Estimating the width of the outlying times
DataWidth = zeros(size(DataBreakIndices));
for m=1:1:size(DataBreakIndices)
    if DataBreakIndices(m) == 1
        DataWidth(m) = round((MCS.TimeStamp(DataBreakIndices(m))).*3600 ./ ...
                             (MCS.ProfilesPerHistogram(DataBreakIndices(m))./LaserRepRate));
    else
        DataWidth(m) = round((MCS.TimeStamp(DataBreakIndices(m))-MCS.TimeStamp(DataBreakIndices(m)-1)).*3600 ./ ...
                             (MCS.ProfilesPerHistogram(DataBreakIndices(m)-1)./LaserRepRate));
    end
    
end

% Making sure to only include complete scans
DataWidth = round(DataWidth./size(A,1)).*size(A,1);
              
%% Converting the names structure to a cell array
CellArray  = struct2cell(MCS);
FieldNames = fieldnames(MCS);

%% Filling data arrays (if needed)
if isempty(B)
    NewCell = CellArray;
else
    % Pre-allocating data cell array
    NewCell = cell(size(CellArray));
    % Looping over cells to fill them with empty data
    for n=1:1:size(CellArray,1)
        %%%%%%%% Initializing data array %%%%%%%%
        NewCell{n,1} = [];
        % Initializing the array location
        StartIndex    = 1;
        EndOfDayBreak = 0;
        for m=1:1:size(B)
            EndIndex = DataBreakIndices(m) - 1;
            %%%%%%%% Determining what to put into the time gaps %%%%%%%%
            if strcmp(FieldNames{n,1},'TimeStamp')
                if EndIndex ~= 0
                    ToAdd = linspace(CellArray{n,1}(EndIndex), ...
                                     CellArray{n,1}(EndIndex+1), ...
                                     DataWidth(m).*size(A,1) + 2)';
                else
                    ToAdd = linspace(0, ...
                                     CellArray{n,1}(EndIndex+1), ...
                                     DataWidth(m).*size(A,1) + 2)';
                end
                ToAdd = ToAdd(2:end-1);
            elseif strcmp(FieldNames{n,1},'Channel')
                ToAdd = repmat(A,DataWidth(m),1);
            elseif strcmp(FieldNames{n,1},'Data')
                ToAdd = repmat(A.*nan,DataWidth(m),size(CellArray{n,1},2));
            else
                ToAdd = repmat(A.*nan,DataWidth(m),1);
            end
            %%%%%%%% Inserting data into time gaps %%%%%%%%
            if DataBreakIndices(m) == 1
                %%%%%%%% There is a break at the start of the day %%%%%%%%
                % Determining if there are multiple breaks in the day,
                % otherwise, just make the end index the last measurement
                if size(DataBreakIndices,1)>1
                    EndIndex = DataBreakIndices(m+1)-1;
                else
                    EndIndex = size(CellArray{n,1},1);
                end
                % Padding the new cell array
%                 NewCell{n,1} = [NewCell{n,1};
%                                 ToAdd;
%                                 CellArray{n,1}(StartIndex:EndIndex,:)];
%                 StartIndex = EndIndex+1;
                NewCell{n,1} = [NewCell{n,1};
                                ToAdd];
%                 StartIndex = EndIndex+1;
            elseif DataBreakIndices(m) == size(CellArray{n,1},1) -1
                %%%%%%%% There is a break at the end of the day %%%%%%%%
                NewCell{n,1} = [NewCell{n,1};
                                CellArray{n,1}(StartIndex:EndIndex,:);
                                ToAdd];
                EndOfDayBreak = 1;
            else
                %%%%%%%% There is a break in the middle of the day %%%%%%%%
                % Insert meaningless time stamps in the middle of the day
                NewCell{n,1} = [NewCell{n,1};
                                CellArray{n,1}(StartIndex:EndIndex,:);
                                ToAdd];
                StartIndex = EndIndex+1;
            end
            %%%%%%%% FIlling end of day if needed %%%%%%%%
            if m == size(B,1) && EndOfDayBreak == 0
                NewCell{n,1} = [NewCell{n,1};
                                CellArray{n,1}(StartIndex:end,:)];
            end
        end
    end
end
%% Converting back to a named structure
MCS = cell2struct(NewCell,FieldNames);
end
