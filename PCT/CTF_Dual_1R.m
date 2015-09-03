function [sino] = CTF_Dual_1R(IR,out,sigma)
% CTF phase retrieval duality method for measured intensity data from one 
% distance only
%       
% In:   IR      intensity measurements
%       out     output from forward model PCT_forward
%       sigma   propotionality constant: -delta/beta
%
% Out: 
%       sino    retrieved mu sinogram
%
% By Rasmus Dalgas Kongskov, 28/05/2014, DTU

% initialiez variables
[Nt,Nth] = size(IR);

% proportionality constant for duality method (sigma)
if nargin<3
    sigma = -out.index{1}(1)/out.index{1}(2);
end

% CTF model
P = 2*sigma*sin(out.w)-2*cos(out.w);
bf = 1/sqrt(Nt)*fft(IR-1);

% mask to avoid near-zero crossing amplification
Z  = abs(P)>0;

sino = sqrt(Nt)*ifft(Z.*(bf./P));

end