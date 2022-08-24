


function [Raw,Data1D,Scan,Availible] = IdentifyNeededInfo(Data,Cal,As,Chan,Const)
%
%
%
%
%% Assuming the data will be found
Availible = true;
%%
for m=1:1:length(As)
    % Loading count profiles
    try
        Raw.(As{m}).TimeStamp = Data.Lidar.Interp.([As{m},Chan{m}]).TimeStamp.*60.*60;
        Raw.(As{m}).Range     = (0:size(Data.Lidar.Interp.([As{m},Chan{m}]).Data,2)-1)'...
                                 .*mean(Data.Lidar.Interp.([As{m},Chan{m}]).RangeResolution.*Const.C*1e-9/2,'omitnan');
        Raw.(As{m}).Value = Data.Lidar.Interp.([As{m},Chan{m}]).Data';
        % Loading wavelength info
        Data1D.Wavelength.(As{m}).TimeStamp = Data.TimeSeries.Laser.(As{m}).TimeStamp.*60.*60;
        Data1D.Wavelength.(As{m}).Value     = Data.TimeSeries.Laser.(As{m}).WavelengthActual./1e9;
        % Interpolating wavelength to handle wavemeter dropouts
        Data1D.Wavelength.(As{m}).Value = InterpolateWavelength(Data1D.Wavelength.(As{m}).TimeStamp,...
                                                                Data1D.Wavelength.(As{m}).Value);
        % Loading MCS info
        Data1D.MCS.(As{m}) = Data.Lidar.Interp.([As{m},Chan{m}]);
        Data1D.MCS.(As{m}) = rmfield(Data1D.MCS.(As{m}),'Data');
        Data1D.MCS.(As{m}).TimeStamp = Data1D.MCS.(As{m}).TimeStamp.*60.*60;
        % Loading wavelength scan info
        Scan.(As{m}) = Cal.ScanData.([As{m},Chan{m}]);
    catch
        Raw = []; Data1D = []; Scan = [];
        Availible = false;
        break
    end
end
end

function [WLReturn] = InterpolateWavelength(Time,WL)
%
%
%
%% Removing bad data (i.e. data equal to nan)
Time2 = Time;
Time2(isnan(WL)) = [];
WL(isnan(WL)) = [];
%% Performing interpolation of data and excluding all extrapolation
WLReturn = interp1(Time2,WL,Time,'linear',nan);
end