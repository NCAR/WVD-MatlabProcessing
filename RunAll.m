
close all; clear; clc;

%% Options to define what systems and when to analyze MPD data
Dates2Process   = {'20211116','20211206'};
Systems2Process = {'mpd_01'};

%% Options to define what to process
ProcessHousekeeping = false;
ProcessRetrievals   = true;

%% Processing all required data
DateNum = datenum(Dates2Process,'yyyymmdd');
TStart = tic;
for n=DateNum(1):1:DateNum(2)
    for m=1:1:length(Systems2Process)
        DateString = datestr(n,'yyyymmdd');
        fprintf(['Processing ',upper(erase(Systems2Process{m},'_')),': ',DateString,' (',datestr(now,'HH:MM:SS'),')\n'])
        try
            [~,~,~,~,~] = RunLoader(DateString,Systems2Process{m},'Skinny',ProcessHousekeeping,ProcessRetrievals);
        catch
            fprintf('************Processing Failed************\n')
        end
    end
end
TEnd = toc(TStart);
fprintf('Elapsed time is %0.2f minutes.\n',TEnd/60)
