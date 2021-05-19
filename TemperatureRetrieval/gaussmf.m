function [f] = gaussmf(x,c,sigma)
f = exp(-(x-c).^2./(2.*sigma.^2));
end