% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August 20th, 2020

function [Thinned] = PoissonThinLidarData(RawData,BinInfo,Options)
%
% Inputs: RawData: Structure containing raw lidar data from a given channel
%         BinInfo: A structure containing the user defined binning settings
%         Options: User defined options for any particular processing 
%
% Output: Thinned: Strucutre containining original data, thinned data, and 
%                  background subtracted data
%
%% Converting structure to cell array for iterating
[Cell,FieldNames,~] = RecursiveStruct2Cell(RawData);
%% Looping over range and time to bin data
for m=1:1:size(Cell,1)
    % Pulling out channel data from structure and thinning
    CData = RawData.(FieldNames{m,1});
    P = PythonThinning(CData);
    % Background subtracting the Poisson thinned data
    for n=1:1:size(P,1)
        BSubCell{n,1}{m,1} = BGSubtractLidarData(P{n,1},CData,BinInfo,Options);
        % Modifying the field names cell array
        if m == 1
            FN{n,1} = ['PoissThined',num2str(n)];
            FN{n,2} = FieldNames;
        end
    end
end
%% Convert the count cell array back to a structure
Thinned = RecursiveCell2Struct(BSubCell, FN);
end

function [P] = Thinning(CData)
%
% Inputs: CData: Structure containing all lidar data from a given channel
%
% Outputs: P:    Cell array containing 2 Poisson thinned arrays
%
%%
P{1,1} = arrayfun(@(x) binornd(x,0.5),CData.Counts);
P{2,1} = CData.Counts - P{1,1};
end

function [P] = PythonThinning(CData)
%
% Inputs: CData: Structure containing all lidar data from a given channel
%
% Outputs: P:    Cell array containing 2 Poisson thinned arrays
%
%%
%% Finding where photon counts should be masked (can't deal with negatives)
Mask = (CData.Counts<0) | isnan(CData.Counts);
CData.Counts(Mask) = 1;      % Removing negative values
%% Poisson thinning
P{1,1} = double(pyrunfile('PoissThin.py','F',x=int64(CData.Counts)));
P{1,1}(Mask) = nan;
P{2,1} = CData.Counts - P{1,1};
end

