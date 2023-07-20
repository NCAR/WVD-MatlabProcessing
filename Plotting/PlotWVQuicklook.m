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
%% Checking for weird figure settings
Options = OverwriteColorbar(Options,Op.Date,Op.System);
%% Setting up needed variables
Sc = [1,1,2560,1416];
Date = datestr(datenum(Op.Date,'yyyymmdd'), 'dd mmm yyyy');
%% Loading colormap 
RBColormap = importdata('CM_NCAR.mat');
%% Checking for python data (to add Cramer Rao lower bound info to plot)
[AddCRLB,Py] = RecursivelyCheckIsField(WV, 'Python');
if AddCRLB; AddCRLB = not(isequal(size(Py.TimeStamp),[1,1])); end
%% Formatting the figure as a whole
FigNum = FindCurrentFigure + 1;
figure(FigNum)
set(gcf,'Position',[Sc(3)/20 Sc(4)/10 Sc(3)/1.5 Sc(4)/1.5],'PaperUnits', 'points', 'PaperPosition', [0 0 1280 800]);
%% Plotting Relative Backscatter
if AddCRLB; Subs = 5; Range = 1:2; else; Subs = 4; Range = 1:2; end
Ax(1) = subplot(Subs,1,Range);
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
Ax(2) = subplot(Subs,1,Range);
pcolor(WV.TimeStamp./60./60,WV.Range./1e3,WV.Smoothed2);
colormap(gca,'jet')
% Formatting the subplot
FormatAxis(Date, Options.WV, 'Water Vapor (g/m^3)',Subs,Range);
%% Plotting Cramer Rao Lower Bound online wavelength
if AddCRLB
   % Mask to remove data without valid relative backscatter in the column
    A = double(all(isnan(WV.RB.Value)));
    Mask = interp1(WV.RB.TimeStamp,A,Py.TimeStamp) > 0;
    if all(size(Mask) == size(Py.ActualWavelength))
        Py.ActualWavelength(Mask) = nan;
    else
        Py.ActualWavelength = nan.*Mask;
    end
    if all(size(Mask) == size(Py.OptimalWavelength))
        Py.OptimalWavelength(Mask) = nan;
    else
        Py.OptimalWavelength = nan.*Mask;
    end
    % Plotting the wavelengths of interst
    Ax(3) = subplot(5,1,5); hold on;
    plot(Py.TimeStamp./60./60,Py.ActualWavelength,'k',Py.TimeStamp./60./60,Py.OptimalWavelength,'r','linewidth',2);
    % Labeling the plot
    xlabel('Time (UTC)'); ylabel('Wavelength (nm)');
    legend({'Observed','CRLB'},'location','south','orientation','horizontal','AutoUpdate','off');
    % Setting the plot bounds
    xlim([0,24]); set(gca,'xtick',0:2:24); YL = ylim; MinYL = 0.02;
    if (YL(2)-YL(1)) < MinYL; ylim(([-MinYL,MinYL]./2)+mean(YL)); end
end
%% Formatting figure
FormatFigures(gcf);
%% Making it so the figure doesn't resize when saving
set(gcf, 'PaperPositionMode', 'auto');
%% Setting the axis location to compensate for width of colorbars
if AddCRLB
    for m=1:1:length(Ax)
        Pos(m,:) = get(Ax(m),'position');
    end
    Width = min(Pos(:,3));
    for m=1:1:length(Ax)
        set(Ax(m),'position',[Pos(m,1),Pos(m,2),Width,Pos(m,4)],'LineWidth',2);
    end
end
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
    title({[Date,' MPD Data',Op.ExtraName]});
end
%% Formatting axes
set(gca,'xtick',0:2:24,'ytick',0:2:Op.MaxRange/1e3)
%% Setting axis bounds
xlim([0,24]); ylim([Op.MinRange,Op.MaxRange]./1e3);
%% Formatting the colorbar
caxis(Op.CAxis); CB = colorbar('EastOutside');
ylabel(CB,Name)
end

function Options = OverwriteColorbar(Options, Date, System)
%
%
%
%
%%
if any(strcmp(System,{'mpd_02','mpd_03','mpd_04'}))            &&  ...
    datenum(Date,'yyyymmdd') >= datenum('20220520','yyyymmdd') &&  ...
    datenum(Date,'yyyymmdd') <= datenum('20220901','yyyymmdd')
    switch System
        case 'mpd_02'; Loc = 'NCU';
        case 'mpd_03'; Loc = 'Yilan';
        case 'mpd_04'; Loc = 'Hsinchu';
        otherwise; Loc = '';
    end
    % Precip
    Options.WV.CAxis = [0,25];
    Options.RB.ExtraName = [' (PRECIP -- ',Loc,')'];
    Options.WV.ExtraName = '';
elseif any(strcmp(System,{'mpd_02','mpd_03'}))                 &&  ...
    datenum(Date,'yyyymmdd') >= datenum('20230712','yyyymmdd') &&  ...
    datenum(Date,'yyyymmdd') <= datenum('20230905','yyyymmdd')
    % M2HATS
    Loc = 'Tonopah';
    Options.WV.CAxis = [0,8];
    Options.RB.ExtraName = [' (M2HATS -- ',Loc,')'];
    Options.WV.ExtraName = '';
else
    Options.RB.ExtraName = '';
    Options.WV.ExtraName = '';
end
end
