% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August 20th, 2020

function [InterpData] = BinLidarData(RawData,Afterpulse,Op)
%
%
%
%
%
%
%% Parsing out options
Times          = Op.TimeGridLidar;
DefaultOptions = Op.Default;
%% Defining processing variables
ResCheck  = {'RangeResolution';'NBins'};       % Variables to define resolution
ToAverage = {ResCheck{:},'RTime'}'; %#ok<CCAT> % Variables to average over 
ToSum     = {'Data';'ProfilesPerHistogram'};   % Variables to sum over
%% Applying the afterpulse correction to the raw data
if Op.Afterpulse
for m = fields(RawData)'
    % Interpolating the afterpulse contour to current range grid if needed
    if size(Afterpulse.(m{1}).Afterpulse,1) == size(RawData.(m{1}).Data,2)
        AP = Afterpulse.(m{1}).Afterpulse'      .* ... % Rate     [Hz ]
             RawData.(m{1}).ProfilesPerHistogram.*...  % Shots    [   ]
             (RawData.(m{1}).RangeResolution.*1e-9);   % BinWidth [sec]
    else
        fprintf(['!!!!!!Failed to apply afterpulse correction: ',m{1},'!!!!!!\n'])
        AP = zeros(size(RawData.(m{1}).Data));
    end
   RawData.(m{1}).Data = RawData.(m{1}).Data - AP;
end
end
%% Converting count structure to cell array
FieldNames = fieldnames(RawData);
RawData    = struct2cell(RawData);
%% Checking if all data has the same range and time resolution
% Checking for unique elements of ResCheck in all Data types
UniqueRes   = cellfun(@(R) cellfun(@(X) unique(R.(X)),ResCheck,'Uni',false),RawData,'Uni',false);
% Checking if there is more than 1 resolution in any data or data type
Sizes       = cellfun(@(B) cellfun(@(A) size(A,1),B),UniqueRes,'Uni',false);
% Setting the needs interpolation flag 
NeedsInterp = size(unique(vertcat(Sizes{:})),1) > 1;
%% Pushing data to know grid by looping over data types
InterpData = cell(size(RawData)); % Pre-allocating cell array
Start = (1:size(Times)-1)';       % Defining counter to get the start times 
for m=1:1:size(RawData,1)
    % Determining the time stamps as the middle of the data period
    InterpData{m,1}.TimeStamp = mean([Times(Start),Times(Start+1)],2,'omitnan');
    % Performing interpolation if needed
    if NeedsInterp
        RawData{m,1}=InterpDataToDefaultGrid(RawData{m,1},DefaultOptions);
    end
    % Looping over interpolation time frames
    for n=1:1:size(Times)-1
        % Finding all data in the desired time frame
        Elements = find(RawData{m,1}.TimeStamp >= Times(n) & ...
                        RawData{m,1}.TimeStamp <  Times(n+1));
        % Combining data on the identified time grid
        for p=1:1:size(ToAverage,1)
            InterpData{m,1}.(ToAverage{p})(n,:) = mean(RawData{m,1}.(ToAverage{p})(Elements,:),'omitnan');
        end
        % Summing variables to be combined
        for p=1:1:size(ToSum,1)
            InterpData{m,1}.(ToSum{p})(n,:) = sum(RawData{m,1}.(ToSum{p})(Elements,:),'omitnan');
        end         
    end
    % Removing data that has not yet come in time
    % Summing variables to be combined
    for p=1:1:size(ToSum,1)
        InterpData{m,1}.(ToSum{p})(InterpData{m,1}.(ToSum{p}) == 0) = nan;
    end
end
%% Converting back to a structure
InterpData = cell2struct(InterpData,FieldNames);
end

function [Raw] = InterpDataToDefaultGrid(Raw,Defaults)
%
%
%
%
%%
% Loop over all data

% Determine if data needs to be binned in range to be lower res than default

% Interpolate data

end