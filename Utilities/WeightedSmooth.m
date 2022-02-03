% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created February, 2022

function [Smoothed] = WeightedSmooth(Data,Options)
%
% Inputs: Data:      An array of anything to be smoothed
%         Options:   A structure containing the user defined options
%                    defining smoothing parameters
%
% Outputs: Smoothed: An array of anything with smoothing applied
%
%% Making the smoothing kernal based on user defined options
SKern  = MakeSmoothingKernal(Options);
%% Running a weighted smoothing algorithm
% Determining where weights should be ignored
B = ~isnan(Data);   
% Removing bad data
Data(~B) = 0;        
% Running the smoothing filter knowing the bad data is in there and then 
% dividing by the fraction of good data
Smoothed = filter2(SKern,Data,'same')./filter2(SKern,B,'same');
Smoothed(~B) = nan;
end

function [Smoothing] = MakeSmoothingKernal(Options)
%
% Inputs: Options:    A structure containing the user defined options
%                     defining smoothing parameters
%
% Outputs: Smoothing: A 2d gaussian smoothing kernal array
%
%% Doubling width for smoothing window & making the FWHM half of it
[X,Y] = meshgrid(linspace(-1,1,2.*Options.SmoothTime/Options.BinTime),...
                 linspace(-1,1,2.*Options.SmoothRange/Options.BinRange)); 
Smoothing = Gaussian2D(X,Y,0,0,0.5,0.5);
end

function [Value] = Gaussian2D(x,y,xo,yo,Sx,Sy)
%
% Inputs: x:  X location of the smoothing kernal (normalized units)
%         y:  Y location of the smoothing kernal (normalized units)
%         xo: X location of the center point of the smoothing kernal
%         yo: Y location of the center point of the smoothing kernal
%         Sx: Standard deviation of the smoothing kernal in the x direction
%         Sy: Standard deviation of the smoothing kernal in the y direction
%
% Outputs: Value: A 2d gaussian smoothing kernal array
%
%% Calculating a normalized 2d Gaussian smoothing kernal
Value = exp(-((((x-xo).^2)/(2.*Sx.*Sx))+(((y-yo).^2)/(2.*Sy.*Sy))));
Value = Value./sum(sum(Value));
end