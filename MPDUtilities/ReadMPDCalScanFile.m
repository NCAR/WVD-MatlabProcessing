% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created January, 2020

function [S] = ReadMPDCalScanFile(FieldNames,File)
%
% Inputs: FieldNames: Cell array of field names of the found photon count
%                     MCS data used to define what scan fields to load
%         File:       String defining full path of calibration scan file 
%
% Outputs: S:         Structure containing all of the needed scan data
%
%% Defining file element names 
ToLoad = {'_scan_wavelength','_Transmission'};
SaveAs = {'Wavelength',      'Transmission'};
%% Loading data
for F = FieldNames(:)'
    for m=1:1:length(ToLoad)
        S.(F{1,1}).(SaveAs{m}) = ReadVariable(File,[F{1,1},ToLoad{m}]);
    end
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