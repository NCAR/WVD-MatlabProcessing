


function [Smoothing] = MakeSmoothingKernal(Options)
%
%
%
%
%% Doubling width for smoothing window & making the FWHM half of it
[X,Y] = meshgrid(linspace(-1,1,2.*Options.SmoothTime/Options.BinTime),...
                 linspace(-1,1,2.*Options.SmoothRange/Options.BinRange)); 
Smoothing = Gaussian2D(X,Y,0,0,0.5,0.5);
end

function [Value] = Gaussian2D(x,y,xo,yo,Sx,Sy)
%
%
%
%
%%
Value = exp(-((((x-xo).^2)/(2.*Sx.*Sx))+(((y-yo).^2)/(2.*Sy.*Sy))));
Value = Value./sum(sum(Value));
end