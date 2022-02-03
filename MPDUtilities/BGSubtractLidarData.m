% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created February 3rd, 2023

function [BGS] = BGSubtractLidarData(Counts,CData,BinInfo,Op)
%
% Inputs: Counts:  An array containing the counts to be background
%                  subtracted. This can be equal to the counts in the CData
%                  structure or Poisson thinned
%         CData:   Structure containing all of the lidar data from a
%                  particular channel (original data with no thinning)
%         BinInfo: A structure containing the user defined binning settings
%         Op:      User defined options for any particular processing 
%
% Output: BGS:     A strucutre containins original lidar data plus the
%                  background subtracted lidar data
%
%% Determining the background index
A = floor(Op.BackgroundInd/BinInfo.BinNum(contains(BinInfo.BinDir,'Range')));
%%  Background subtracting but only for data to max altitude
BGS.Background = mean(Counts(end-A:end,:));
BGS.Counts     = Counts(CData.Range<Op.MaxRange,:) - BGS.Background;
BGS.Raw        = Counts(CData.Range<Op.MaxRange,:);
BGS.Range      = CData.Range(CData.Range<Op.MaxRange);
BGS.TimeStamp  = CData.TimeStamp;
end