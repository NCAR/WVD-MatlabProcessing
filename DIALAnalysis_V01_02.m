% Modified by Stillwell in Jan-Feb 2018
        % Clear out all unused (commented out) lines and fat 
        % Pushing sub-tasks into subfunctions for readability
        % Changed variable names to push data into structures/cell arrays

        % This function is to replicate Scott's code with nothing added
        % except reading new data format (use with new Labview)
% Modified by Stillwell in Oct-Nov 2018    
        % Compartmentalized all retrievals
        % Renewed ability to run multiple wavelengths per day
        % Reorganizing code into subfunctions more rigidly 
        % Commenting code more thoroughly for other users
        
function [Uptime] = DIALAnalysis_V01_02(JSondeData, Options, Paths)
%
% Inputs: JSondeData: A structure containing all of the loaded calibration
%                     data from the JSonde files
%         Options:    A structure containing all of the user defined
%                     processing options
%         Paths:      A structure containing all of the relevant file path
%                     information 
%
% Outputs: none
%
%%
tic;
fprintf(['Processing: ',Options.System,' data from 20',Paths.Date,'\n'])

%% Reading Jsond file information
[Capabilities,Iterations,HardwareMap,Map] = ReadFakeJSondFile(Paths);

DataTypes = {'Etalon*.nc';'LL*.nc';'MCS*.nc';'Pow*.nc';
             'WS*.nc';'UPS*.nc';'HKeep*.nc';'Humidity*.nc'};

%% Importing all netcdf data files from the selected date
% Loading data
[Counts,PulseInfoNew,PulseInfo.Uptime] = ReadRawNetCDFData(DataTypes,HardwareMap,Paths);
Uptime = PulseInfo.Uptime;
% Determining pulse info
PulseInfo.BinWidth    = round((double(nanmean(PulseInfoNew.Data.RangeResolution{1,1}))*1e-9*3e8/2)*10)/10;
PulseInfo.DataTimeRaw = double(PulseInfoNew.TimeStamp.LidarData{1,1})./24 + ...
                        day(datetime(['20',Paths.Date],'inputformat','yyyyMMdd'),'dayofyear');
PulseInfo.DeltaRIndex = ceil(75/PulseInfo.BinWidth); % this is the cumlative sum photons gate spacing
PulseInfo.DeltaR      = PulseInfo.DeltaRIndex*PulseInfo.BinWidth*100; % delta r in cm
time_per_column       = nanmedian(PulseInfoNew.Data.ProfilesPerHistogram{1,1})./7e3;
PulseInfo.Profiles2AverageWV = 2*round(((Options.ave_time.wv*60/time_per_column)+1)/2); % 7kHz, 10k accum data rate is ~1.4s
PulseInfo.Profiles2AverageRB = 2*round(((Options.ave_time.rb*60/time_per_column)+1)/2); % 7kHz, 10k accum data rate is ~1.4s
PulseInfo.DataTimeInt = downsample(movmean(PulseInfo.DataTimeRaw,PulseInfo.Profiles2AverageRB),PulseInfo.Profiles2AverageRB); 
PulseInfo.DataTimeInt = PulseInfo.DataTimeInt(1:floor(size(PulseInfo.DataTimeRaw,1)/PulseInfo.Profiles2AverageRB));
clear time_per_column DataTypes
%% Defining arrays used to smooth the data
AverageRange   = [1;round(1500/PulseInfo.BinWidth);round(2500/PulseInfo.BinWidth)];
SpatialAverage = [150/PulseInfo.BinWidth; 300/PulseInfo.BinWidth; 600/PulseInfo.BinWidth];

%% Recast weather station data
if Options.flag.WS==1
  SurfaceWeather          = PulseInfoNew.WeatherStation;
  SurfaceWeather.Pressure = SurfaceWeather.Pressure./1013.249977;  % pressure in atm
  % Calculate surface parameters from weatherstation data
  [~,SurfaceWeather.NumberDensity] = ConvertWeatherStationValues(SurfaceWeather.RelativeHumidity,SurfaceWeather.Temperature);
end

%% Initial data preparation
fprintf('Processing Data\n')
% grid data in time to final array size
PulseInfo.DataTime = (floor(min(PulseInfo.DataTimeRaw)):1/24/60*(Options.ave_time.gr):(floor(min(PulseInfo.DataTimeRaw))+1))';
% remove the time lag from cumsum
PulseInfo.DataTimeShifted = PulseInfo.DataTime -(1/24/60*((Options.ave_time.wv-1)/2));
%Calculating average wavelength
for m=1:1:size(Counts.Raw,1)
    PulseInfoNew.Laser.WavelengthActual{m,1} = PulseInfoNew.Laser.WavelengthActual{m,1} - JSondeData.WavemeterOffset;
    PulseInfoNew.Laser.WavelengthDesired{m,1} = PulseInfoNew.Laser.WavelengthDesired{m,1} - JSondeData.WavemeterOffset;
end
clear m

%% Checking for multiple wavelengths
for m=1:1:size(PulseInfoNew.Laser.WavelengthDesired,1)
   PulseInfo.LambdaNearest{m,1} = round(PulseInfoNew.Laser.WavelengthDesired{m,1},4);
   PulseInfo.LambdaUnique{m,1}  = unique(round(PulseInfoNew.Laser.WavelengthDesired{m,1}(...
                                              isnan(PulseInfoNew.Laser.WavelengthDesired{m,1}) == 0),4));
end
  
%% Range vector in meters
Altitude.RangeAcquired = single(0:PulseInfo.BinWidth:(size(Counts.Raw{1,1},2)-1)*PulseInfo.BinWidth);
if PulseInfo.DeltaRIndex ~= 1
    Altitude.RangeOriginal = downsample(movmean(Altitude.RangeAcquired,PulseInfo.DeltaRIndex),PulseInfo.DeltaRIndex,1);
else
    Altitude.RangeOriginal = Altitude.RangeAcquired;
end
Altitude.RangeShift    = (PulseInfo.DeltaRIndex-1)/2*PulseInfo.BinWidth + JSondeData.RangeCorrection; %
Altitude.RangeActual   = Altitude.RangeOriginal+Altitude.RangeShift; % actual range points of data

%% Processing photon counts
i = size(Counts.Raw{1,1}, 1);         % Number of time bins
j = size(Altitude.RangeAcquired, 2);  % Number of altitude bins

% Looping over channels and performing operations on photon counting data
for m=1:1:size(Counts.Raw,1)
    fprintf(['Processing ',Map.Channels{m},' count arrays\n'])
    % Applying saturation correction
    if Options.flag.pileup == 1
        Counts.Parsed{m,1} = CorrectPileUp(single(Counts.Raw{m,1}),PulseInfoNew.Data.ProfilesPerHistogram{m,1}, ...
                                               PulseInfoNew.Data.RangeResolution{m,1},JSondeData.DeadTime);
    else
        Counts.Parsed{m,1} = single(Counts.Raw{m,1});
    end
    % Integrating photon counts in time and space
    Counts.Integrated{m,1} = Bin(Counts.Parsed{m,1},i,j,PulseInfo.Profiles2AverageRB,PulseInfo.DeltaRIndex);
    % select last ~1200 meters to measure background
    Counts.Background1D{m,1} = mean(Counts.Integrated{m,1}(:,end-round(1200/PulseInfo.BinWidth./PulseInfo.DeltaRIndex):end),2);
    % Background subtracting the parsed counts
    Counts.BackgroundSubtracted{m,1} = (bsxfun(@minus, Counts.Integrated{m,1}, Counts.Background1D{m,1}));
    % smooth RB for 1 minute and set spatial average 
    Counts.RelativeBackscatter{m,1} = movmean(Counts.BackgroundSubtracted{m,1},2,1,'omitnan');
    % Range correcting relative backscatter
    Temp = PulseInfoNew.Data.RangeResolution{m,1}.* PulseInfoNew.Data.ProfilesPerHistogram{m,1}.*(1-JSondeData.SwitchRatio);
    Temp = downsample(movmean(Temp,PulseInfo.Profiles2AverageRB),PulseInfo.Profiles2AverageRB);  
    Temp = Temp(1:floor(i/PulseInfo.Profiles2AverageRB));
    Counts.RelativeBackscatter{m,1} = bsxfun(@times,  Counts.RelativeBackscatter{m,1},Altitude.RangeOriginal.^2./Temp./PulseInfo.Profiles2AverageRB./PulseInfo.DeltaRIndex);
    % Converting background to count rate (counts/sec
    Counts.Background1D{m,1} = Counts.Background1D{m,1}./(Temp.*1e-9); 
    % Integrating photon counts in time and space
%    size(Counts.Integrated{m,1})
%    size(Counts.Parsed{m,1})
%    Counts.Integrated{m,1}(isnan(Counts.Parsed{m,1})) = nan;
    % Regular averaging
    Counts.CountRate{m,1}    = Counts.BackgroundSubtracted{m,1}./PulseInfo.Profiles2AverageRB./PulseInfo.DeltaRIndex;
    % Grid to regular gate spacing in time (???????make recursive???????)
    Counts.Background1D{m,1}        = interp1(PulseInfo.DataTimeInt, Counts.Background1D{m,1}        ,PulseInfo.DataTime,Options.InterpMethod,Options.Extrapolation);
    Counts.CountRate{m,1}           = interp1(PulseInfo.DataTimeInt, Counts.CountRate{m,1}           ,PulseInfo.DataTime,Options.InterpMethod,Options.Extrapolation);
    Counts.BackgroundSubtracted{m,1}= interp1(PulseInfo.DataTimeInt, Counts.BackgroundSubtracted{m,1},PulseInfo.DataTime,Options.InterpMethod,Options.Extrapolation);
    Counts.RelativeBackscatter{m,1} = interp1(PulseInfo.DataTimeInt, Counts.RelativeBackscatter{m,1} ,PulseInfo.DataTime,Options.InterpMethod,Options.Extrapolation);
    % Removing any negative counts
    Counts.CountRate{m,1}(real(Counts.CountRate{m,1}) <= 0) = 0;
    Counts.RelativeBackscatter{m,1}(real(Counts.RelativeBackscatter{m,1}) <= 0) = 0;
end
clear m Temp i j


%% Gradient filter for the WV data
if Options.flag.gradient_filter == 1
    [FX,~] = gradient(Counts.CountRate{Map.Offline,1});
    % Determining the allowable limit of the gradient
    Limit = (PulseInfoNew.Data.RangeResolution{Map.Offline,1}./250.* PulseInfoNew.Data.ProfilesPerHistogram{Map.Offline,1}./14000).*1000;
    Limit = downsample(movmean(Limit,PulseInfo.Profiles2AverageRB),PulseInfo.Profiles2AverageRB);  
    Limit = Limit(1:floor(size(PulseInfoNew.Data.RangeResolution{Map.Offline,1},1)/PulseInfo.Profiles2AverageRB));
    Limit = interp1(PulseInfo.DataTimeInt,Limit,PulseInfo.DataTime,Options.InterpMethod,Options.Extrapolation);
    Counts.CountRate{Map.Offline,1}(FX<-Limit | FX> Limit) = nan; % remove falling (leading) edge of clouds
    clear FX Limit
end

%% Recursively interpolate data from collected to averaged grid
if Options.flag.WS == 1
    SurfaceWeather = RecursivelyInterpolateStructure(SurfaceWeather,PulseInfo.DataTimeRaw,PulseInfo.DataTime,Options.InterpMethod,Options.Extrapolation);
else
    SurfaceWeather.Pressure    = [];
    SurfaceWeather.Temperature = [];
end
Temp = PulseInfo.DataTimeRaw;
PulseInfo = RecursivelyInterpolateStructure(PulseInfo,double(PulseInfo.DataTimeRaw),double(PulseInfo.DataTime),Options.InterpMethod,Options.Extrapolation);
PulseInfo.DataTimeRaw = Temp;
clear Temp
%% Performing data retrievals
DataProducts = []; % Pre-allocating retrieval info 
for m=1:1:Iterations
    fprintf('   Retrieval Iteration %0.0f\n',m)
    %%%%%%%%%%% Performing the water vapor retrievals %%%%%%%%%%%
    if Capabilities.WVDIAL == 1
        [Counts,DataProducts] = RetrievalsH2O(Altitude,Counts,DataProducts,   ...
                                              JSondeData,Map,Options,Paths,   ...
                                              PulseInfo,PulseInfoNew,         ...
                                              SpatialAverage, SurfaceWeather, ...
                                              AverageRange);
    end
    %%%%%%%%%%% Performing the HSRL retrievals %%%%%%%%%%%
    if Capabilities.HSRL == 1 || Capabilities.O2HSRL == 1
        % RetrievalsHSRL(Altitude,Capabilities,Counts,DataProducts,Map,Options,Paths,PulseInfo,SurfaceWeather);
    end
    %%%%%%%%%%% Performing perterbative temperature retrievals %%%%%%%%%%%
    if Capabilities.O2DIAL == 1
        
    end
end
  
%% Finding plotting information
PulseInfo.DataTimeDateNumFormat = datenum(2000+str2double(Paths.Date(1:2)),1,0)+ double(PulseInfo.DataTime);

%% Decimate data in time to final array size
if Options.flag.decimate == 1
    decimate_time  = 1;%Options.ave_time.wv/Options.ave_time.gr; %ave_time.wv/ave_time.gr;
    if size(Altitude.RangeOriginal,2) == 280
        decimate_range = 1; % keep native gate spacing
    else
        decimate_range = 2; % push all data to lower resolution
    end
    % Decimating needed count profiles
    for m=1:1:size(Counts.RelativeBackscatter,1)
        % average RB data before decimating
        Counts.RelativeBackscatter{m,1} = movmean(Counts.RelativeBackscatter{m,1},decimate_time*2,1,'omitnan');
        % Decimating
        Counts.RelativeBackscatter{m,1} = Counts.RelativeBackscatter{m,1}(1:decimate_time:end, 1:decimate_range:end);
        Counts.Background1D{m,1}        = Counts.Background1D{m,1}(       1:decimate_time:end, 1:decimate_range:end);
        Counts.CountRate{m,1}           = Counts.CountRate{m,1}(          1:decimate_time:end, 1:decimate_range:end);
    end
    % Decimating altitude to the final grid spacing
    Altitude.RangeOriginal = Altitude.RangeOriginal(1:decimate_range:end);
    % Decimate data products in space and time
    DataProducts = RecursivelyDecimateStructure(DataProducts,decimate_time ,size(DataProducts.OpticalDepth,1), ...
                                                             decimate_range,size(DataProducts.OpticalDepth,2));    
    % Decimating pulseinfo and surface weather time series data in time
    PulseInfo = RecursivelyDecimateStructure(PulseInfo,decimate_time,size(PulseInfo.DataTime,1),[],[]);        
    if Options.flag.WS ==1
        SurfaceWeather = RecursivelyDecimateStructure(SurfaceWeather,decimate_time,size(SurfaceWeather.AbsoluteHumidity,1),[],[]);        
    end
end

%% Plot Data
Plotting = PlotData(Altitude,Counts,DataProducts,Map,Options,Paths,PulseInfo,PulseInfoNew,SurfaceWeather); %#ok<NASGU>

%% Save Data
fprintf('Saving Data\n')
if Options.flag.save_data == 1
    cd(Paths.SaveData)
    Paths.FileName = ['ProcessedMPDData_MPD0',Options.System(6),'_20',num2str(Paths.Date),'.mat'];
    save(Paths.FileName,'Altitude','Counts','DataProducts','Options','Paths','Plotting','PulseInfo','PulseInfoNew','SurfaceWeather','-v7.3')
    cd(Paths.Code)
end
if Options.flag.save_netCDF == 1  % save the data as an nc file
    WriteNetCDFData (PulseInfo.LambdaMedian{2,1} ,PulseInfo.LambdaMedian{1,1},DataProducts.N_avg,DataProducts.N_Error,                  ...
                     Counts.CountRate{1,1},Counts.CountRate{2,1}, ...
                     P,Counts.RelativeBackscatter{1,1},Altitude.RangeOriginal,T,PulseInfo.DataTimeDateNumFormat)
end

%% Cleaning the workspace variables that are unneeded
clear decimate_range decimate_time AverageRange SpatialAverage m
toc
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Integrated] = Bin(CountArray,TimeBins,AltBins,Prof2Avg,Heights2Avg)
%
%
%
%
%% Integrating in time
Temp  = cumsum(CountArray,1,'omitnan')-[zeros(Prof2Avg,AltBins);cumsum(CountArray(1:TimeBins-Prof2Avg,:),1,'omitnan')];  %rolling average of rows or time
%% Integrating in range
Integrated  = cumsum(Temp,2,'omitnan')-[zeros(TimeBins,Heights2Avg), cumsum(Temp(:,1:AltBins-Heights2Avg),2,'omitnan')]; % rolling sum of collumns or range
%% Only returning the completely interpolated count arrays
Integrated = Integrated(Prof2Avg:Prof2Avg:end,Heights2Avg:Heights2Avg:end);
    
end

function [Surf_AH,Surf_N] = ConvertWeatherStationValues(Surf_RH,Surf_T) 
%
% Inputs: Surf_RH:  Surface relative humidity time series   [Units: %]
%         Surf_T:   Surface temperature time series         [Units: C]
%                   
% Outputs: Surf_AH: Surface absolute humiditiy time series  [Units: g/kg]
%          Surf_N:  Surface number density time series      [Units: 1/cm^3]
%
%% Defining universal constants
[Constants,Conversions] = DefineConstants;

%% Defining conversion constants for Clausius Clapeyron Equation
a0 = 6.107799961;
a1 = 4.436518521E-1;
a2 = 1.428945805E-2;
a3 = 2.650648471E-4;
a4 = 3.031240396E-6;
a5 = 2.034080948E-8;
a6 = 6.136820929E-11;
%% Calculating vapor pressure of water
e=((a0+Surf_T.*(a1+Surf_T.*(a2+Surf_T.*(a3+Surf_T.*(a4+Surf_T.*(a5+Surf_T.*a6))))))./1); %vapor pressure in hPa
%% Convert RH to number density and absolute humidity
Surf_N  = (Surf_RH.*e./(Constants.R.*Conversions.CelciusToKelvin(Surf_T))).*Constants.N_A*1e-6; %cm^-3
Surf_AH = Surf_N.*Constants.m.*1e9;   % g/m^3  
end

function [CorrCounts] = CorrectPileUp(Counts,MCSAccum,BinDuration,t_d)
%
% Inputs: Counts:      
%         MCS: 
%         t_d
%                      
% Outputs: CorrCounts: 
%
%%
% MCSC gives counts accumulated for set bin duration so convert to count rate  C/s.
% divide by bin time in sec (e.g., 500ns) and # of acumulations (e.g., 10000)
% e.g., 10 accumlated counts is 2000 C/s
CorrFactor = 1./(1-(t_d.*(Counts./(BinDuration.*1E-9.*MCSAccum))));
CorrCounts = Counts.*CorrFactor;
end

function WriteNetCDFData (lambda,lambda_off,N_avg,N_error,Offline_Temp_Spatial_Avg,Online_Temp_Spatial_Avg,P,RB,range,T,time_new)
%
% Inputs: lambda:
%         lambda_off:
%         N_avg:
%         N_error:
%         Offline_Temp_Spatial_Avg:
%         Online_Temp_Spatial_Avg
%         P:
%         RB:
%         range:
%         T
%         time_new:
%
%% Convert NaN fill values (and Inf) to -1 flag
time_new(isnan(time_new)==1) = -1;  
N_avg(isnan(N_avg)==1)       = -1;
N_error(isnan(N_error)==1)   = -1;
N_avg(isinf(N_avg)==1)       = -1;
N_error(isinf(N_error)==1)   = -1;
Offline_Temp_Spatial_Avg(isinf(Offline_Temp_Spatial_Avg)==1) = -1;
Online_Temp_Spatial_Avg(isinf(Online_Temp_Spatial_Avg)==1)   = -1;
time_unix = (time_new-datenum(1970,1,1))*86400; % convert to unix time
   
%% Creating the NetCDF file and dimensions
% Define the file
cdf_name = strcat('wv_dial.', datestr(date, 'yymmdd'), Paths.FolderType);
ncid     = netcdf.create([cdf_name '.nc'],'CLOBBER');
% define the dimensions and variables
dimid1 = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
dimid2 = netcdf.defDim(ncid, 'range', length(range));
dimid3 = netcdf.defDim(ncid, 'lambda', length(lambda));
    
%% Defining the variables
myvarID1 = netcdf.defVar(ncid,'time','double',dimid1);
  netcdf.putAtt(ncid, myvarID1, 'units', 'days since January 0, 0000')
myvarID2 = netcdf.defVar(ncid,'range','float',dimid2);
  netcdf.putAtt(ncid, myvarID2, 'units', 'meters')
myvarID3 = netcdf.defVar(ncid,'N_avg','float',[dimid2 dimid1]);
  netcdf.putAtt(ncid, myvarID3, 'long_name', 'water_vapor_number_density')
  netcdf.putAtt(ncid, myvarID3, 'units', 'molecules/cm^3')
  netcdf.putAtt(ncid, myvarID3, 'FillValue', '-1')
myvarID4 = netcdf.defVar(ncid,'N_error','float',[dimid2 dimid1]);
  netcdf.putAtt(ncid, myvarID4, 'long_name', 'water_vapor_number_density_error')
  netcdf.putAtt(ncid, myvarID4, 'units', 'molecules/cm^3')
  netcdf.putAtt(ncid, myvarID4, 'FillValue', '-1')
myvarID5 = netcdf.defVar(ncid,'P','float', dimid2);
  netcdf.putAtt(ncid, myvarID5, 'long_name', 'pressure')
  netcdf.putAtt(ncid, myvarID5, 'units', 'atm')
myvarID6 = netcdf.defVar(ncid,'T','float', dimid2);
  netcdf.putAtt(ncid, myvarID6, 'long_name', 'temperature')
  netcdf.putAtt(ncid, myvarID6, 'units', 'degK')
myvarID7 = netcdf.defVar(ncid,'RB','float',[dimid2 dimid1]);
  netcdf.putAtt(ncid, myvarID7, 'long_name', 'relative_backscatter')
  netcdf.putAtt(ncid, myvarID7, 'units', 'arbitrary_units')
  netcdf.putAtt(ncid, myvarID7, 'FillValue', '-1')
myvarID8 = netcdf.defVar(ncid,'time_unix','double',dimid1);
  netcdf.putAtt(ncid, myvarID8, 'long_name', 'unix time')
  netcdf.putAtt(ncid, myvarID8, 'units', 'seconds since 00:00:00 UTC, January 1, 1970')
  netcdf.putAtt(ncid, myvarID8, 'FillValue', '-1')
myvarID9 = netcdf.defVar(ncid,'Offline_Temp_Spatial_Avg','float',[dimid2 dimid1]);
  netcdf.putAtt(ncid, myvarID9, 'long_name', 'offline counts')
  netcdf.putAtt(ncid, myvarID9, 'units', 'counts')
  netcdf.putAtt(ncid, myvarID9, 'FillValue', '-1')
myvarID10 = netcdf.defVar(ncid,'Online_Temp_Spatial_Avg','float',[dimid2 dimid1]);
  netcdf.putAtt(ncid, myvarID10, 'long_name', 'online counts')
  netcdf.putAtt(ncid, myvarID10, 'units', 'counts')
  netcdf.putAtt(ncid, myvarID10, 'FillValue', '-1')
myvarID11 = netcdf.defVar(ncid,'lambda','float', dimid3);
  netcdf.putAtt(ncid, myvarID11, 'long_name', 'online wavelength')
  netcdf.putAtt(ncid, myvarID11, 'units', 'nm')
  netcdf.putAtt(ncid, myvarID11, 'FillValue', '-1')
myvarID12 = netcdf.defVar(ncid,'lambda_off','float', dimid3);
  netcdf.putAtt(ncid, myvarID11, 'long_name', 'offline wavelength')
  netcdf.putAtt(ncid, myvarID11, 'units', 'nm')
  netcdf.putAtt(ncid, myvarID11, 'FillValue', '-1')
netcdf.endDef(ncid)

%% Save the variables to the file
netcdf.putVar(ncid,myvarID1,0,length(time_new),time_new);
netcdf.putVar(ncid,myvarID2,range);
netcdf.putVar(ncid,myvarID3,[0,0],size(N_avg'),N_avg');
netcdf.putVar(ncid,myvarID4,[0,0],size(N_error'),N_error');
netcdf.putVar(ncid,myvarID5,P);
netcdf.putVar(ncid,myvarID6,T);
netcdf.putVar(ncid,myvarID7,[0,0],size(RB'),RB');
netcdf.putVar(ncid,myvarID8,time_unix);
netcdf.putVar(ncid,myvarID9,[0,0],size(Offline_Temp_Spatial_Avg'),Offline_Temp_Spatial_Avg');
netcdf.putVar(ncid,myvarID10,[0,0],size(Online_Temp_Spatial_Avg'),Online_Temp_Spatial_Avg');
netcdf.putVar(ncid,myvarID11,lambda);
netcdf.putVar(ncid,myvarID12,lambda_off);
netcdf.close(ncid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
