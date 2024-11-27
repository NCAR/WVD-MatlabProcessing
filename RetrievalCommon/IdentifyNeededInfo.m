


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
        % Checking if the channel is low or not
        if contains(As{m},'Low')
            SaveAs = As{m}(1:end-3);
        else
            SaveAs = As{m};
        end
        Raw.(SaveAs).TimeStamp = Data.Lidar.Interp.([As{m},Chan{m}]).TimeStamp.*60.*60;
        Raw.(SaveAs).Range     = (0:size(Data.Lidar.Interp.([As{m},Chan{m}]).Data,2)-1)'...
                                 .*mean(Data.Lidar.Interp.([As{m},Chan{m}]).RangeResolution.*Const.C*1e-9/2,'omitnan');
        Raw.(SaveAs).Value = Data.Lidar.Interp.([As{m},Chan{m}]).Data';
        % Loading wavelength info
        Data1D.Wavelength.(SaveAs).TimeStamp = Data.TimeSeries.Laser.(SaveAs).TimeStamp.*60.*60;
        Data1D.Wavelength.(SaveAs).Value     = Data.TimeSeries.Laser.(SaveAs).WavelengthActual./1e9;
        % Interpolating wavelength to handle wavemeter dropouts
        Data1D.Wavelength.(SaveAs).Value = InterpolateWavelength(Data1D.Wavelength.(SaveAs).TimeStamp,...
                                                                 Data1D.Wavelength.(SaveAs).Value);
        % Loading wavelength scan info
        Scan.(SaveAs) = Cal.ScanData.([SaveAs,Chan{m}]);
        % Loading MCS info
        Data1D.MCS.(SaveAs) = Data.Lidar.Interp.([As{m},Chan{m}]);
        Data1D.MCS.(SaveAs) = rmfield(Data1D.MCS.(SaveAs),'Data');
        Data1D.MCS.(SaveAs).TimeStamp = Data1D.MCS.(SaveAs).TimeStamp.*60.*60;
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