% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August 20th, 2020

function [Thinned] = PoissonThinLidarData(RawData,BinInfo,Options)
%
%
%
%
%
%
%% Converting structure to cell array for iterating
[Cell,FieldNames,~] = RecursiveStruct2Cell(RawData);
%% Determining the background index
A = floor(Options.BackgroundInd/BinInfo.BinNum(contains(BinInfo.BinDir,'Range')));
%% Looping over range and time to bin data
for m=1:1:size(Cell,1)
    % Pulling out channel data from structure
    ChannelData = RawData.(FieldNames{m,1});
    % Performing the thinning on a single profile
    P{1,1} = arrayfun(@(x) binornd(x,0.5),ChannelData.Counts);
    P{2,1} = ChannelData.Counts - P{1,1};
    % Background subtracting
    for n=1:1:size(P,1)
        % Calculating background
        BSubCell{n,1}{m,1}.Background = mean(P{n,1}(end-A:end,:));
        % Background subtracting but only for data to max altitude
        BSubCell{n,1}{m,1}.Counts     = P{n,1}(ChannelData.Range<Options.MaxRange,:) - BSubCell{n,1}{m,1}.Background;
        BSubCell{n,1}{m,1}.Raw        = P{n,1}(ChannelData.Range<Options.MaxRange,:);
        BSubCell{n,1}{m,1}.Range      = ChannelData.Range(ChannelData.Range<Options.MaxRange);
        BSubCell{n,1}{m,1}.TimeStamp  = ChannelData.TimeStamp;
        % Modifying the field names cell array
        if m == 1
            FN{n,1} = ['PoissThined',num2str(n)];
            FN{n,2} = FieldNames;
        end
    end
    % Cleaning workspace
    clear ChannelData P
end
%% Convert the count cell array back to a structure
Thinned = RecursiveCell2Struct(BSubCell, FN);
end

