% Written By: Robert Stillwell
% Written For: NCAR

function TopPlot(X,Y,Color,Style,Width)
hold on;
h = plot(X,Y,'Color',Color,'LineStyle',Style,'LineWidth',Width);
set(h,'marker','none')
uistack(h,'top')
end