% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created March, 2022

function FigNum = PlotWVQuicklook(WV,Options,Op)
%
% Inputs: WV:
%         Options: Options structure for just WV
%         Op: Full options structure
%
% Outputs: none
%
%% Setting up needed variables
Sc   = get(0,'ScreenSize');
Date = datestr(datenum(Op.Date,'yyyymmdd'), 'dd mmm yyyy');
%% Loading colormap 
RBColormap = importdata('NCAR_C_Map.mat');
%% Checking for python data (to add Cramer Rao lower bound info to plot)
[AddCRLB,Py] = RecursivelyCheckIsField(WV, 'Python');
if AddCRLB; AddCRLB = not(isequal(size(Py.TimeStamp),[1,1])); end
%% Formatting the figure as a whole
FigNum = FindCurrentFigure + 1;
figure(FigNum)
set(gcf,'Position',[Sc(3)/20 Sc(4)/10 Sc(3)/1.5 Sc(4)/1.5],'PaperUnits', 'points', 'PaperPosition', [0 0 1280 800]);
%% Plotting Relative Backscatter
if AddCRLB; Subs = 5; Range = 1:2; else; Subs = 4; Range = 1:2; end
subplot(Subs,1,Range);
pcolor(WV.RB.TimeStamp./60./60,WV.RB.Range./1e3,real(log10(WV.RB.Value*5)));
colormap(gca,RBColormap)
% Formatting the subplot
cb = FormatAxis(Date, Options.RB, 'Relative Backscatter (C/ns km^2)',Subs,Range);
% Reformatting colorbar ticks to make user know it is logrithmic
set(cb,'ytick',Options.RB.CAxis(1):1:Options.RB.CAxis(2))
A = get(cb,'yticklabel'); C = get(cb,'ytick');
for m=1:1:size(A,1); B{m,1} = ['10^{',A{m},'}']; end      %#ok<AGROW>
set(cb,'ytick',C,'yticklabel',B)
%% Plotting Water Vapor
% Applying the water vapor data mask
WV.Smoothed2(WV.Mask) = nan;
% Plotting
if AddCRLB; Subs = 5; Range = 3:4; else; Subs = 4; Range = 3:4; end
subplot(Subs,1,Range);
pcolor(WV.TimeStamp./60./60,WV.Range./1e3,WV.Smoothed2);
colormap(gca,'jet')
% Formatting the subplot
FormatAxis(Date, Options.WV, 'Water Vapor (g/m^3)',Subs,Range);
%% Plotting Cramer Rao Lower Bound online wavelength
if AddCRLB
    % Getting handle of WV contour
    Ax1 = gca;
    % Plotting the wavelengths of interst
    subplot(5,1,5); hold on;
    plot(Py.TimeStamp./60./60,Py.ActualWavelength,'k',Py.TimeStamp./60./60,Py.OptimalWavelength,'r','linewidth',2);
    % Labeling the plot
    xlabel('Time (UTC)'); ylabel('Wavelength (nm)');
    legend({'Observed','CRLB'},'location','north','orientation','horizontal','AutoUpdate','off');
    % Setting the plot bounds
    xlim([0,24]); set(gca,'xtick',0:2:24); YL = ylim; MinYL = 0.015;
    if (YL(2)-YL(1)) < MinYL; ylim(([-MinYL,MinYL]./2)+mean(YL)); end
    % Getting handle of CRLB plot
    Ax2 = gca;
end
%% Formatting figure
FormatFigures;
%% Setting the axis location to compensate for width of colorbars
if AddCRLB
    Pos = get(Ax1,'position'); Pos2 = get(Ax2,'position');
    set(Ax2,'position',[Pos2(1),Pos2(2),Pos(3),Pos2(4)],'LineWidth',2);
end
%% Making it so the figure doesn't resize when saving
set(gcf, 'PaperPositionMode', 'auto');
end

function CB = FormatAxis(Date, Op, Name,Subplots,Range)
%
%
%
%
%
%% Formatting the pcolor type
shading interp
%% Labeling plot
if any(ismember(Range,Subplots))
    xlabel('Time (UTC)');
end
ylabel('Height (km, AGL)');
if any(ismember(Range,1))
    title({[Date,' MPD Data']});
end
%% Formatting axes
set(gca,'xtick',0:2:24,'ytick',0:2:Op.MaxRange/1e3)
%% Setting axis bounds
xlim([0,24]); ylim([Op.MinRange,Op.MaxRange]./1e3);
%% Formatting the colorbar
caxis(Op.CAxis); CB = colorbar('EastOutside');
ylabel(CB,Name)
end

