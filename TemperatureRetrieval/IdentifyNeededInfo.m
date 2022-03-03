


function [Raw,Data1D,Scan,Availible] = IdentifyNeededInfo(Data,Cal)
%
%
%
%
%% Assuming the data will be found
Availible = true;
%%
As   = {'O2Online';'O2Offline'};
Chan = {'Comb';    'Comb'};
for m=1:1:length(As)
    % Loading count profiles
    try
        Raw.(As{m}).TimeStamp = Data.Lidar.Interp.([As{m},Chan{m}]).TimeStamp.*60.*60;
        Raw.(As{m}).Range     = (0:size(Data.Lidar.Interp.([As{m},Chan{m}]).Data,2)-1)'...
                                 .*mean(Data.Lidar.Interp.([As{m},Chan{m}]).RangeResolution.*299792458*1e-9/2,'omitnan');
        Raw.(As{m}).Value = Data.Lidar.Interp.([As{m},Chan{m}]).Data';
        % Loading wavelength info
        Data1D.Wavelength.(As{m}).TimeStamp = Data.TimeSeries.Laser.(As{m}).TimeStamp.*60.*60;
        Data1D.Wavelength.(As{m}).Value = Data.TimeSeries.Laser.(As{m}).WavelengthActual./1e9;
        % Laoding wavelength scan info
        Scan.(As{m}) = Cal.ScanData.([As{m},Chan{m}]);
    catch
        Raw = []; Data1D = []; Scan = [];
        Availible = false;
        break
    end
end
end