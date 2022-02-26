


function [Raw,Data1D,Availible] = IdentifyNeededInfo(Data)
%
%
%
%
%% Assuming the data will be found
Availible = true;
%%
As = {'O2Online';'O2Offline'};
for m=1:1:length(As)
    % Loading count profiles
    try
        Raw.(As{m}).TimeStamp = Data.Lidar.Interp.([As{m},'Comb']).TimeStamp.*60.*60;
        Raw.(As{m}).Range     = (0:size(Data.Lidar.Interp.([As{m},'Comb']).Data,2)-1)'...
                                 .*mean(Data.Lidar.Interp.([As{m},'Comb']).RangeResolution.*299792458*1e-9/2,'omitnan');
        Raw.(As{m}).Value = Data.Lidar.Interp.([As{m},'Comb']).Data';
        % Loading wavelength info
        Data1D.Wavelength.(As{m}).TimeStamp = Data.TimeSeries.Laser.(As{m}).TimeStamp.*60.*60;
        Data1D.Wavelength.(As{m}).Value = Data.TimeSeries.Laser.(As{m}).WavelengthActual./1e9;
    catch
        Raw = []; Data1D = [];
        Availible = false;
        break
    end
end
end