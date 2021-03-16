
close all; clear; clc;

%% Options to define what systems and when to analyze MPD data
% Dates2Process   = {'20200829','20200910'};
Dates2Process   = {'20200829','20200829'};
% Dates2Process   = {'20200829','20210131'};
% Dates2Process   = {'20200829','20200829'};

Systems2Process = {'mpd_01','mpd_02','mpd_03','mpd_04','mpd_05'};
% Systems2Process = {'mpd_01'};

%% Processing all required data
DateNum = datenum(Dates2Process,'yyyymmdd');
tic
for n=DateNum(1):1:DateNum(2)
    for m=1:1:length(Systems2Process)
        DateString = datestr(n,'yyyymmdd');
        fprintf(['Processing ',upper(erase(Systems2Process{m},'_')),': ',DateString,' (',datestr(now,'HH:MM:SS'),')\n'])
        try
            [~,~,~,~,~] = RunLoader(DateString,Systems2Process{m},'None');
        catch
            fprintf('************Processing Failed************\n')
        end
    end
end
toc