


function FormatFigures()
%
%
%
%%
Figures = sort(double(findall(0,'type','figure')));
for m=1:1:length(Figures)
    set(0,'CurrentFigure',Figures(m));
    AddSolidOutlines(gcf)
end
end


function AddSolidOutlines(CurrentFigure)
%
% Inputs: CurrentFigure: Output of the gcf function
%
%% Plotting constants
Pts       = 1000;
GridColor = [1,1,1]*0.75;
GridWidth = 0.75;
BoxWidth  = 2;
FontSize  = 16;
%% Finding all the subplots
subs=sort(unique(findall(findall(CurrentFigure),'Type','axes','Tag','')));
%% Formatting subplots
for n=1:1:size(subs,1)
    % Setting the active axis
    axes(subs(n))  %#ok<LAXES>
    set(subs(n),'FontSize',FontSize);
    % Determining the locations of the needed lines
    XLims = xlim; YLims = ylim;
    XTickLoc = get(gca,'xtick'); 
    YTickLoc = get(gca,'ytick'); 
    % Turning the grid off because we are re-creating it
    grid off;
    % Plotting vertical lines
    for m=1:1:size(XTickLoc,2)
        TopPlot(XTickLoc(m)*ones(1,Pts),linspace(YLims(1),YLims(2),Pts),GridColor,'--',GridWidth)
    end
    % Plotting horizontal lines
    for m=1:1:size(YTickLoc,2)
        TopPlot(linspace(XLims(1),XLims(2),Pts),YTickLoc(m)*ones(1,Pts),GridColor,'--',GridWidth)
    end
    % Plotting outside box lines
    for m=1:1:2
        % Horizontal box lines
        TopPlot(linspace(XLims(1),XLims(2),Pts),YLims(m)*ones(1,Pts),'k','-',BoxWidth)
        % Vertical box lines
        TopPlot(XLims(m)*ones(1,Pts),linspace(YLims(1),YLims(2),Pts),'k','-',BoxWidth)
    end
    % Making sure axis doesn't resize
    axis tight
end
end
