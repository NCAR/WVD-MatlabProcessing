% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created July 31st, 2020

function PlotTSData(RawData,Data,FigNum)
%
% Inputs: RawData: A data structure containing loaded data
%         Data:    A structure containing organized and interpolated data
%         FigNum:  A number of the current figure being plotted 
%
%% Plotting constants
FigLoc   = [10,92,1002,746];
%% Checking inputs
if nargin == 2
    FigNum = FindCurrentFigure;
end
%% Parsing out all levels of data structure
[Data2Plot{1,1},Data2Plot{1,2}] = FindDataStructures(RawData);
[Data2Plot{2,1},Data2Plot{2,2}] = FindDataStructures(Data);
%% Defining the line type to plot
Type = {'k-';'g.'};   % First is raw and second is interpolated
%% Looping over raw and interpolated data and plotting
% FiguresMade = [];
for n=1:1:size(Data2Plot,1)
    % Looping over elements of the data structure
    for m=1:1:size(Data2Plot{n,1},1)
        % Pulling out data from data structure
        TimeStamp  = Data2Plot{n,1}{m,1}.TimeStamp;
        ToPlot     = rmfield(Data2Plot{n,1}{m,1},'TimeStamp');
        FieldNames = fieldnames(ToPlot);
        ToPlot     = struct2cell(ToPlot);
        SubPlots   = size(ToPlot,1);
        % Plotting all data sources as a subplot
        figure(m + FigNum)
%         Handle = figure(m + FigNum);
%         FiguresMade(end+1) = Handle;
        set(gcf,'position',FigLoc,'Color',[1 1 1])
        for p=1:1:SubPlots
            subplot(SubPlots,1,p); hold on;
            plot(TimeStamp,ToPlot{p,1},Type{n,1},'linewidth',2);
            ylabel(FieldNames{p,1})
            xlim([min(TimeStamp),max(TimeStamp)])
            if p == 1
               title([Data2Plot{n,2}{m,1},' Data']) 
            elseif p == SubPlots
                xlabel('Time Stamp [UTC]'); 
            end
            if mod(p,2) == 0
               set(gca,'YAxisLocation','right');
            end
        end
    end
end
end

