

function [Data1D,Data2D,Found] = LoadPythonData(FileName,Op)
%
%
%
%
%% Checking if python file exists...if not escape
Found = false;
[~, FName, FExt] = fileparts(FileName);
if isfile(FileName)
    [~, FName, FExt] = fileparts(FileName);
    CWLogging(['    Loading: ',[FName,FExt],'\n'],Op,'Main')
else
    Data1D = []; Data2D = [];
    CWLogging(['    ** Error: ',[FName,FExt],' can''t be found.\n'],Op,'Main')
    return
end
%% Loading NCIP Data
% Map from file names (Types) to code names (As)
Types    = {'Temperature_Model';'Pressure_Model'};
As       = {'Temperature';'Pressure'};
% Map from file variables (SubTypes) to code variables (SubAs)
ForTypes = {   'time'  ; 'range';   ''  };
AftTypes = {     ''    ;   ''   ;   ''  };
SubAs    = {'TimeStamp';'Range' ;'Value'};
% Loading data
Data2D.NCIP = LoadData(Types,ForTypes,AftTypes,As,SubAs,FileName);
%% Loading Surface Weather Station Data
% Map from file names (Types) to code names (As)
Types    = {'Surface_Temperature';'Surface_Pressure'};
As       = {'Temperature';'Pressure'};
% Map from file variables (SubTypes) to code variables (SubAs)
ForTypes = {   'time'  ;   ''  };
AftTypes = {     ''    ;   ''  };
SubAs    = {'TimeStamp';'Value'};
% Loading data
Data1D = LoadData(Types,ForTypes,AftTypes,As,SubAs,FileName);
%% Loading retrieved MPD data
% Map from file names (Types) to code names (As)
Types    = {'Backscatter_Ratio';'Absolute_Humidity';'Aerosol_Backscatter_Coefficient'};
As       = {'BackRatio';'Humidity';'BSCoefficient'};
% Map from file variables (SubTypes) to code variables (SubAs)
ForTypes = {   'time'  ; 'range';   ''  ;     ''    ;  ''   };
AftTypes = {     ''    ;   ''   ;   ''  ;'_variance';'_mask'};
SubAs    = {'TimeStamp';'Range' ;'Value';'Variance' ;'Mask' };
% Loading data
Data2D.MPD = LoadData(Types,ForTypes,AftTypes,As,SubAs,FileName);
%% Reading range offsets (and removing if it can be found)
Res = ReadAttribute(FileName,'range','range_resolution');
Off = ReadAttribute(FileName,'range','bin_offset');
if not(isnan(Res*Off))
    Data2D.MPD = RecursiveOverwriteField(Data2D.MPD,'Range',[],-Res*Off);
end
%% 
Found = true;
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

% This function tries to read an attribute of interest. If the type is Double
% it will typecast the variable, otherwise the type is maintained from the
% netcdf file.
function [A] = ReadAttribute(Filename,FileVar,Att,VariableType)
%
% Inputs: Filename:     String containing the desired file name
%         FileVar:      NetCDF variable name to check
%         Att:          Attribute of the FileVar to load
%         VariableType: Data type of variable to read
%
% Outputs: A:           Loaded data
%
%% Checking data inputs
if nargin == 4
    VariableType = 'Double';
end
%% Loading data
A = [];
try
    if strcmp(VariableType,'Double')
        A = double(ncreadatt(Filename,FileVar,Att));
    elseif strcmp(VariableType,'Native')
        A = ncreadatt(Filename,FileVar,Att);
    end
catch
    A = nan;
end
end
