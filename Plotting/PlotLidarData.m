% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August 15th, 2020

function PlotLidarData(Data,FigNum)
%
% Inputs: RawData: A data structure containing loaded data
%         Data:    A structure containing organized and interpolated data
%         FigNum:  A number of the current figure being plotted 
%
%% Constants
C = 299792458;
%% Plotting constants
FigLoc   = [10,92,1105,713];
CAxisLim = [1,7];
%% Checking inputs
if nargin == 1
    FigNum = FindCurrentFigure;
end
%% Parsing out all levels of data structure
[Data2Plot{1,1},Data2Plot{1,2}] = FindDataStructures(Data);
%% Looping over raw and interpolated data and plotting
for n=1:1:size(Data2Plot,1)
    % Looping over elements of the data structure
    for m=1:1:size(Data2Plot{n,1},1)
        % Pulling out data from data structure
        ToPlot     = Data2Plot{n,1}{m,1};
        TimeStamp  = Data2Plot{n,1}{m,1}.TimeStamp;
        ToPlot     = rmfield(ToPlot,'TimeStamp');
        CountData  = Data2Plot{n,1}{m,1}.Data;
        ToPlot     = rmfield(ToPlot,'Data');
        RangeRes   = Data2Plot{n,1}{m,1}.RangeResolution;
        ToPlot     = rmfield(ToPlot,'RangeResolution');
        % Setting up ancillary data to be plotted
        FieldNames = fieldnames(ToPlot);
        ToPlot     = struct2cell(ToPlot);
        SubPlots   = size(ToPlot,1);
        % Creating range and time arrays
        Altitude = (0:1:size(CountData,2)-1).*(mean(RangeRes,'omitnan')./1e9).*(C./2)./1e3;
        % Plotting count data
        figure(m + FigNum);
        set(gcf,'position',FigLoc,'Color',[1 1 1])
        subplot(SubPlots*2,1,1:SubPlots)
        pcolor(TimeStamp, Altitude, real(log10(CountData')));
        shading flat; grid on; box on; caxis(CAxisLim);
        title([Data2Plot{n,2}{m,1},' Photon Counts']); ylabel('Altitiude [km]');
        % Locating the colorbar
        Pos = get(gca,'position');
        A   = colorbar; CBWidth = get(A,'position'); CBWidth = CBWidth(3);
        set(A,'position',[Pos(1) + Pos(3) + CBWidth./3,Pos(2),CBWidth,Pos(4)])
        CBTicks = get(A,'ytick');
        for p=1:1:length(CBTicks)
           CBTickLabel{p} = ['10^{',num2str(CBTicks(p)),'}'];  %#ok<AGROW>
        end
        set(A,'ytick',CBTicks,'yticklabel',CBTickLabel)
        % Plotting all data sources as a subplot
        for p=1:1:SubPlots
            subplot(SubPlots*2,1,SubPlots+p); hold on;
            plot(TimeStamp,ToPlot{p,1},'k','linewidth',2); grid on; box on;
            ylabel(FieldNames{p,1});
            xlim([min(TimeStamp),max(TimeStamp)])
            if p == SubPlots
                xlabel('Time Stamp [UTC]');
            end
            if mod(p,2) == 0
               set(gca,'YAxisLocation','right');
            end
        end
    end
end
end

