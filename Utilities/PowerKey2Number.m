


function [Position] = PowerKey2Number(Keys)
%
%
%
%
%%
KeysPossible = {'WVOnline';'WVOffline';'HSRL';'O2Online';'O2Offline'};
%%
Position = ones(size(Keys)).*-1;
for m=1:1:length(Keys)
    for n=1:1:length(KeysPossible)
        if strcmp(Keys{m},KeysPossible{n})
            Position(m) = n;
            break
        end
    end
end
end

