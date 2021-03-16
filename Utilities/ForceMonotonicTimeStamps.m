

function [NewTime] = ForceMonotonicTimeStamps(OldTime)
%
%
%
%%
% Checking for daily rollovers
NewTime = cumsum([0;diff(OldTime)]<0).*24 + OldTime;
% Checking for equal time stamps
NewTime = cumsum([-1;diff(NewTime)]==0).*1e-9 + NewTime;
end