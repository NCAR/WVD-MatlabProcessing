% Written By: Robert Stillwell
% Written For: National Center for Atmospheric Research
%
function [Return] = DensityFiltering(Array,Points,Allowable)
%
% Inputs: Array:     An array of data to be filtered by nan density. This
%                    array will be filtered assuming each measurement is a 
%                    column.
%         Points:    The number of points used to determine what the nan
%                    density about a particular point is.
%         Allowable: The maximum allowable density of nans about the data
%                    point of interest
%
% Outputs: Mask:     An array of the size of the input Array with a nan
%                    where density of input nans is above the allowable
%                    level
%
%% Finding the starting point
ToSkip = Points./2;
Start  = floor(ToSkip);
End    = floor(ToSkip);
%% Arrays
Mask = zeros(size(Array,1)-Points+1,size(Array,2)-Points+1); % Pre-allocating data
for n=1:1:Points
    for m=1:1:Points
        Mask = Mask + Array(m:end-Points+m,n:end-Points+n);
    end
end
Mask = Mask./Points./Points;
%% Removing points with an unallowably high density of nans
Return = Mask > Allowable;
Return = [false(Start,size(Mask,2)+Start+End); 
         [false(size(Mask,1),Start),Return,false(size(Mask,1),End)];
          false(End,size(Mask,2)+Start+End)];
end