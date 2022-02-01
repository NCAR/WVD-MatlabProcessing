% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August 20th, 2020

function [Binned,Binning] = PreProcessLidarData(RawData,Options)
%
%
%
%
%
%
%% Converting structure to cell array for iterating
[Cell,FieldNames,~] = RecursiveStruct2Cell(RawData);
%% Looping over range and time to bin data
BinDir = {'Range';'TimeStamp'};
for m=1:1:size(Cell,1)
    % Pulling out channel data from structure
    ChannelData = RawData.(FieldNames{m,1});
    % Determining the number of bins to pull together
    BinNum = cellfun(@(s) floor(mean(diff(Options.(s)))/mean(diff(ChannelData.(s)))),BinDir);
    % Saving a count array for modification
    Counts = ChannelData.Value;
    GoodData = isnan(Counts);
    Counts(isnan(Counts)) = 0;
    for n=1:1:size(BinDir,1)
        % Bin counts
        CumCounts = cumsum(Counts);
        CumGoodData = cumsum(GoodData);
        % Downsampling to bin resolution
        CumCounts   = CumCounts(BinNum(n):BinNum(n):end,:);
        CumGoodData = CumGoodData(BinNum(n):BinNum(n):end,:);
        % Removing the cumsum
        Counts    =  (CumCounts - [zeros(1,size(CumCounts,2));CumCounts(1:end-1,:)])';
        GoodData    =  (CumGoodData - [zeros(1,size(CumGoodData,2));CumGoodData(1:end-1,:)])';
        InterpCell{m,1}.(BinDir{n}) = ChannelData.(BinDir{n})(BinNum(n):BinNum(n):end);
    end
    Counts(Counts==0 & GoodData >= 1) = nan;
    InterpCell{m,1}.Counts = Counts;
end
%% Convert the count cell array back to a structure
Binned = RecursiveCell2Struct(InterpCell, FieldNames);
Binning.BinNum = BinNum;
Binning.BinDir = BinDir;
end

