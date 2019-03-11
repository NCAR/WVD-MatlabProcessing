



function [RBSpectrum] = BuildRBSpectrumXY(Const,X,Y)
%
%
%
%
%
%
%
%% Calcuating Unknown Tenti Model constants
c_int     = 1.0;
c_tr      = 3/2;
gamma_int = c_int/(c_tr+c_int);
rlx_int   = 1.5*Const.ViscosityBulk/(Const.Viscosity*gamma_int);
eukenf    = Const.MAir*Const.TConductivity/(Const.Viscosity*Const.Kb*(c_tr+c_int));

%% Running the Tenti model
[~,RBSpectrum]=crbs6(Y,rlx_int,eukenf,c_int,c_tr,X);

%% Normalizing RB spectrum
RBSpectrum = RBSpectrum./trapz(X,RBSpectrum);

end

function [cohsig,sptsig,c] = crbs6(y,rlx_int,eukenf,c_int,c_tr,xi)
% Written by John Smith
% October 21st, 2010
% University of Colorado at Boulder, CIRES
% John.A.Smith@Colorado.EDU
% MATLAB version 7.10.0.59 (R2010a) 64-bit
% Adapted from "Coherent Rayleigh-Brillouin Scattering"
% by Xingguo Pan

% Computes the coherent and spontaneous
% RBS spectrum given the parameters
% in crbs_molecular using the s6 model
% by G. Tenti, 1974

% Called by: crbs_molecular.m

n_xi=numel(xi);

n=6;
a=zeros(n,n);
b=zeros(n,2);

cohsig=zeros(1,n_xi);
sptsig=zeros(1,n_xi);

cpxunit=sqrt(-1);

%y7=1.5*y;
gamma_int=c_int/(c_tr+c_int);
j020=-y;
%j030=1.5*j020;
j100=-gamma_int*y/rlx_int;
j001=j100*c_tr/c_int;
j100001=j100*sqrt(c_tr/c_int);
j110=j100*5/6+j020*2/3;
j011110=j100*sqrt(5/(8*c_int));
j_nu=0.4*(1.5+c_int)+(3+c_int)/(2*rlx_int)+9*eukenf/(16*rlx_int^2);
j_de=-1+(4/15)*eukenf*(1.5+c_int)+(c_int/3)*eukenf/rlx_int;
j_co=-y*(2*gamma_int/3);
j011=j_co*j_nu/j_de;

%coharea=0;
%sptarea=0;

for i=1:n_xi
    z=xi(i)+y*cpxunit;
	w0=w0_func(z);
	w1=-sqrt(pi)+z.*w0;
	w2=z.*w1;
	w3=-0.5*sqrt(pi)+z.*w2;
	w4=z.*w3;
	w5=-3*sqrt(pi)/4+z.*w4;
	w6=z.*w5;
    
    i0000=w0/(sqrt(pi));
	i0100=(z.*w0-sqrt(pi))*sqrt(2/pi);
	i0001=i0100;
	i0010=(2*w2-w0)/(sqrt(6*pi));
	i1000=i0010;
	i0011=(2*w3-3*w1)/(sqrt(5*pi));
	i1100=i0011;
	i0101=2*w2/sqrt(pi);
	i0110=(-w1+2*w3)/sqrt(3*pi);
	i1001=i0110;
	i0111=(-3*w2+2*w4)*sqrt(2/(5*pi));
	i1101=i0111;
	i1111=(13*w2-12*w4+4*w6)/(5*sqrt(pi));
	%i0002=(-w0+2*w2)/sqrt(3*pi);
	%i0200=i0002;
	%i0211=(-w1+8*w3-4*w5)/sqrt(15*pi);
	%i1102=i0211;
	%i0202=2*(w0-2*w2+2*w4)/(3*sqrt(pi));
	%i0210=(w0+4*w2-4*w4)/(3*sqrt(2*pi));
	%i1002=i0210;
	%i0102=(-w1+2*w3)*sqrt(2/(3*pi));
	%i0201=i0102;
	i1010=(5*w0-4*w2+4*w4)/(6*sqrt(pi));
	i1110=(7*w1-8*w3+4*w5)/sqrt(30*pi);
	i1011=i1110;
    
    %a_factor=1;
    
    
    a(:,1)=-j020*[i0000 i0001 i0011 i0010 0 0]+[cpxunit 0 0 0 0 0];
    a(:,2)=-j020*[i0100 i0101 i0111 i0110 0 0]+[0 cpxunit 0 0 0 0];
    a(:,3)=(j020-j110)*[i1100 i1101 i1111 i1110 0 0]+j011110*[0 0 0 0 -i0100 -i0101]+[0 0 -cpxunit 0 0 0];
    a(:,4)=(j020-j100)*[i1000 i1001 i1011 i1010 0 0]+j100001*[0 0 0 0 -i0000 -i0001]+[0 0 0 -cpxunit 0 0];
    a(:,5)=j100001*[i1000 i1001 i1011 i1010 0 0]+(j001-j020)*[0 0 0 0 i0000 i0001]+[0 0 0 0 cpxunit 0];
    a(:,6)=j011110*[i1100 i1101 i1111 i1110 0 0]+(j011-j020)*[0 0 0 0 i0100 i0101]+[0 0 0 0 0 cpxunit];
    
    b(:,1)=-[i0100 i0101 i0111 i0110 0 0];
    b(:,2)=-[i0000 i0001 i0011 i0010 0 0];
    
    c=linsolve(a,b);
    
    cohsig(i)=c(1,1)*conj(c(1,1));
    sptsig(i)=2*real(c(1,2));
%    coharea=coharea+cohsig(i);
%    sptarea=sptarea+sptsig(i);
end
end

function w = w0_func(z,N)
% Based on:
% FADDEEVA   Faddeeva function
%   W = FADDEEVA(Z) is the Faddeeva function, aka the plasma dispersion
%   function, for each element of Z. The Faddeeva function is defined as:
%
%     w(z) = exp(-z^2) * erfc(-j*z)
%
%   where erfc(x) is the complex complementary error function.
%
%   W = FADDEEVA(Z,N) can be used to explicitly specify the number of terms
%   to truncate the expansion (see (13) in [1]). N = 16 is used as default.
%
%   Example:
%       x = linspace(-10,10,1001); [X,Y] = meshgrid(x,x); 
%       W = faddeeva(complex(X,Y)); 
%       figure; 
%       subplot(121); imagesc(x,x,real(W)); axis xy square; caxis([-1 1]); 
%       title('re(faddeeva(z))'); xlabel('re(z)'); ylabel('im(z)'); 
%       subplot(122); imagesc(x,x,imag(W)); axis xy square; caxis([-1 1]);
%       title('im(faddeeva(z))'); xlabel('re(z)'); ylabel('im(z)'); 
%
%   Reference:
%   [1] J.A.C. Weideman, "Computation of the Complex Error Function," SIAM
%       J. Numerical Analysis, pp. 1497-1518, No. 5, Vol. 31, Oct., 1994 
%       Available Online: http://www.jstor.org/stable/2158232

% Called by: crbs6.m, crbs7.m


if nargin<2, N = []; end
if isempty(N), N = 50; end

w = zeros(size(z)); % initialize output

%%%%%
% for purely imaginary-valued inputs, use erf as is if z is real
idx = real(z)==0; %
w(idx) = exp(-z(idx).^2).*erfc(imag(z(idx)));

if all(idx), return; end
idx = ~idx;

%%%%%
% for complex-valued inputs

% make sure all points are in the upper half-plane (positive imag. values)
idx1 = idx & imag(z)<0;
z(idx1) = conj(z(idx1));

M = 2*N;
M2 = 2*M;
k = (-M+1:1:M-1)'; % M2 = no. of sampling points.
L = sqrt(N/sqrt(2)); % Optimal choice of L.

theta = k*pi/M;
t = L*tan(theta/2); % Variables theta and t.
f = exp(-t.^2).*(L^2+t.^2);
f = [0; f]; % Function to be transformed.
a = real(fft(fftshift(f)))/M2; % Coefficients of transform.
a = flipud(a(2:N+1)); % Reorder coefficients.

Z = (L+1i*z(idx))./(L-1i*z(idx));
p = polyval(a,Z); % Polynomial evaluation.
w(idx) = 2*p./(L-1i*z(idx)).^2 + (1/sqrt(pi))./(L-1i*z(idx)); % Evaluate w(z).

% convert the upper half-plane results to the lower half-plane if necesary
w(idx1) = conj(2*exp(-z(idx1).^2) - w(idx1));

w=-w*sqrt(-1)*pi; % converts to 'w0_func' used in crbs_molecular.m
end

