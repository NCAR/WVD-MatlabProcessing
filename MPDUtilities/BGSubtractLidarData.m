% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created February 3rd, 2023
% This function is used to background subtract the input data. 2 types of
% input are possible, either a count array or a strucutre containing count
% information. The former is indicated by CData taking teh structure
% information whil ethe latter is identified by CData being empty. If CData
% is empty, the Count structure will be recursively background subtracted,
% otherwise only a single background subtraction will happen.

function [BGS] = BGSubtractLidarData(Counts,CData,BinInfo,Op)
%
% Inputs: Counts:  Type1 = An array containing the counts to be background
%                  subtracted. This can be equal to the counts in the CData
%                  structure or Poisson thinned for example
%                  Type2 = A strucutre containing multiple photon counting
%                  data channels that need to be recursively addressed
%         CData:   Type1 = Structure containing all of the lidar data from
%                  a particular channel (data with no thinning for example)
%                  Type2 = An empty array indicating that the Counts
%                  strucutre needs to be split and recursively acted upon
%         BinInfo: A structure containing the user defined binning settings
%         Op:      User defined options for any particular processing 
%
% Output: BGS:     A strucutre containins original lidar data plus the
%                  background subtracted lidar data. If CData has been
%                  empty, this may be a multi-layered structure
%
%% Stand alone operation requires recursive call to handle all channels
if isempty(CData)
    % Converting structure to cell array for iterating
    [Cell,FieldNames,~] = RecursiveStruct2Cell(Counts);
    % Looping over all channels found
    for m=1:1:size(Cell,1)
        % Pulling out channel data from structure
        CData = Counts.(FieldNames{m,1});
        % Background subtracting each channel
        BSubCell{m,1} = BGSubtractLidarData(CData.Counts,CData,BinInfo,Op); %#ok<AGROW>
    end
    % Convert the count cell array back to a structure
    BGS = RecursiveCell2Struct(BSubCell, FieldNames);
else
    % Determining the background index
    A = floor(Op.BackgroundInd/BinInfo.BinNum(contains(BinInfo.BinDir,'Range')));
    %  Background subtracting but only for data to max altitude
    BGS.Background = mean(Counts(end-A:end,:));
    BGS.Counts     = Counts(CData.Range<Op.MaxRange,:) - BGS.Background;
    BGS.Raw        = Counts(CData.Range<Op.MaxRange,:);
    BGS.Range      = CData.Range(CData.Range<Op.MaxRange);
    BGS.TimeStamp  = CData.TimeStamp;
end
end