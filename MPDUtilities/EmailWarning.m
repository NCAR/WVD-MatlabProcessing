

function EmailWarning(Labels,Times,Status,GLTO,System)
%
%
%
%
%
%
%% 
Fault = 2;            % Value to look for to see if status is bad ("Fault")
TimeToLookBack = 0.5; % Hours
Warning = [];         % Pre-allocating an empty string
%% Labels of elements of interest to check
RequiredData  = {'MCS Data','LL Data','WS Data'};
RequiredPower = {'WVOnlineAmpPower','WVOfflineAmpPower','O2OnlineAmpPower','O2OfflineAmpPower'};
BooleanChecks = {'UPS OnBattery','Temp: HSRLOven','HVACCurrent'};
%%
StartTime = GLTO -TimeToLookBack;
EndTime   = GLTO;
%% Checking that all required data is availible
for m=RequiredData
    % Checking for missing data or warnign flags
    [Missing,WarnReturn] = CheckForWarning(Labels,m{1},Fault,Times,Status,StartTime,EndTime);
    % Appending warnings to the warning string if needed
    Warning = AppendToWarningString(Warning,Missing,WarnReturn,m{1},' seems to be missing',StartTime,EndTime);
end
%% Checking that all laser powers are above zero
for m=RequiredPower
    [Missing,WarnReturn] = CheckForWarning(Labels,m{1},Fault,Times,Status,StartTime,EndTime);
    % Appending warnings to the warning string if needed
    Warning = AppendToWarningString(Warning,Missing,WarnReturn,m{1},' seems low',StartTime,EndTime);
end
%%
for m=BooleanChecks
    %
    [Missing,WarnReturn] = CheckForWarning(Labels,m{1},Fault,Times,Status,StartTime,EndTime);
    % Appending warnings to the warning string if needed
    Warning = AppendToWarningString(Warning,Missing,WarnReturn,m{1},' flag warning',StartTime,EndTime);
end

%% Printing Warings if any
if ~isempty(Warning)
    % Outputting the warnings for the logs
    fprintf(Warning)
    % Eamiling warnings to user
    SendEmail('stillwel@ucar.edu',['MPD Status Error: ',System],Warning)
    % SendEmail('7203293917@txt.att.net','MPD Status Error',Warning)
end

end

function [Missing,WarningReturn] = CheckForWarning(Labels,Desired,Fault,Times,Status,StartTime,EndTime)
%
%
%
%% Determining where to look for data in array
Col = strcmp(Labels,Desired);
if any(Col)
    Column = find(Col,1,'first');
else
    Column = [];
end
%% Pulling out data of interest
% Checking if element is present
if ~isempty(Column)
    % Extracting the dta from the column of interst
    T = Times(Column,:); S = Status(Column,:);
    % Seeing if all data for the look back period is "Fault" (i.e. = 2)
    Missing       = false;
    WarningReturn = all(S(T>StartTime & T<EndTime) == Fault);
else
    Missing       = true;
    WarningReturn = true;
end
end

function [Warning] = AppendToWarningString(Warning,Missing,WarningReturn,Element,SpecificWarn,StartTime,EndTime)
%
%
%
%
%
%%
if Missing
    Warning = strcat(Warning,sprintf([Element,' cant be found for %0.02f - %0.02f UTC.'],StartTime,EndTime),'\n   ');
else
    if WarningReturn
        Warning = strcat(Warning,sprintf([Element,SpecificWarn,' for %0.02f - %0.02f UTC.'],StartTime,EndTime),'\n   ');
    end
end
end

function SendEmail(To,Subject,Body)
%
% Inputs: To:      Recipient of the email
%         Subject: Suject of the email
%         Body:    Main text of the email
%
%%
try
    setpref('Internet','SMTP_Server','localhost');
    warning('off','all')
    NewBody = [strrep(Body,'\n',10)];
    sendmail(To,Subject,[NewBody]) %#ok<NBRAK2>
    warning('on','all')
catch
    fprintf('****************Email failed ot send****************\n')
end

end