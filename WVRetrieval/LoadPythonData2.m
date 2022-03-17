

function [Data1D] = LoadPythonData2(FileName)
%
%
%
%
%% Loading NCIP Data
% Map from file names (Types) to code names (As)
Types    = {'WVOnline_'};
As       = {'Online'};
% Map from file variables (SubTypes) to code variables (SubAs)
ForTypes = {   'time'  ;        ''         ;        ''        };
AftTypes = {     ''    ;'optimalwavelength';'Wavelength_mean' };
SubAs    = {'TimeStamp';'OptimalWavelength';'ActualWavelength'};
% Loading data
Data1D = LoadData(Types,ForTypes,AftTypes,As,SubAs,FileName);

end

function [Data] = LoadData(Types,ForTypes,AftTypes,As,SubAs,FileName)
%
%
%
%% Pre-allocating cell array for loading
Data = cell(size(Types,1),1);
%% Loading data
for m=1:1:size(Types,1)
   % Determining the names of the variables to be read
   Type = cellfun(@(a,b) [a,Types{m,1},b],ForTypes,AftTypes,'Uni',false);
   % Loading each variable
   SubData = cell(length(Type),1);
   for n=1:1:length(Type)
       if strcmp(ForTypes{n},'time')
           SubData{n,1} = ReadVariable(FileName,'time');
       elseif strcmp(ForTypes{n},'range')
           SubData{n,1} = ReadVariable(FileName,'range');
       else
           SubData{n,1} = ReadVariable(FileName,Type{n,1});
       end
   end
   % Converting sub cell array into a structure
   Data{m,1} = cell2struct(SubData,SubAs);
end
%% Converting data from a cell array to structure
Data = cell2struct(Data,As); 
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
    end
catch
    A = nan;
end
end

