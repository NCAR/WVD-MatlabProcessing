% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created February, 2022

function [System] = ReadMPDJSonFiles(File,DDate)
%
% Inputs: File:    Full file path to the desired json file
%         DDate:   String containing desired date to pull file info from
%                  with the format of 'yyyymmdd'
%
% Outputs: System: Structure containing needed MPD information
%
%% Loading Json data
JSon=loadjson(File,'SimplifyCell',1);
%% Defining file information desired
% Structure element in JSon where the data is found
DesiredField = {'Receiver_Scan_File','Wavemeter_offset','Location','Location','Location','Location','Location'};
% Substructure field name where the data is stored
DesiredValue = {'value','value','location','elevation','latitude','longitude','Project'};
% How the user desires saving the fields
ToNameField  = {'ScanFile','WMOffset','Site','SiteAlt','SiteLat','SiteLon','Project'};
% Default values to use
DefaultVal   = {'Unknown',0,'Assumed Boulder',1613,40.0386,-105.2390,{'none'}};
%% Looping and loading data
for m=1:1:length(DesiredField)
    % Extracting structure information to use quickly
    Dates  = extractfield(JSon.(DesiredField{m}),'date');
    Value  = extractfield(JSon.(DesiredField{m}),DesiredValue{m});
    try
       Change = extractfield(JSon.(DesiredField{m}),'change_type'); 
    catch
       Change = extractfield(JSon.(DesiredField{m}),'change_0x20_type'); 
    end    
    % Determining dates of information and how far away those dates are
    Dates = datetime(Dates,'InputFormat','d-MMM-yyyy HH:mm');
    Away  = datenum(DDate,'yyyymmdd') - datenum(Dates);
    % Checking if the data is valid (note that data marked as abrupt must
    % have a time-to-date that is positive, otherwise mark as invalid)
    Away(strcmp(Change,'abrupt') & Away < 0) = nan;
    % Looking for the closest data and saving it
    if isempty(find(abs(Away) == min(abs(Away)),1))
        System.(ToNameField{m}) = DefaultVal{m};
    else
        System.(ToNameField{m}) = Value(abs(Away) == min(abs(Away)));
    end
end
end