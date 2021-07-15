


function [Smoothed] = WeightedSmooth(Data,SmoothKern)
%
%
%
%
%
%% Running a weighted smoothing algorithm
% Determining where weights should be ignored
B = ~isnan(Data);   
% Removing bad data
Data(~B) = 0;        
% Running the smoothing filter knowing the bad data is in there and then 
% dividing by the fraction of good data
Smoothed = filter2(SmoothKern,Data,'same')./ ...
            filter2(SmoothKern,B,'same');
Smoothed(~B) = nan;
end