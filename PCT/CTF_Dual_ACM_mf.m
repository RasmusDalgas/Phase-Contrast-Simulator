    function [A,b,AT] = CTF_Dual_ACM_mf(IR,out,R,sigma)
% Set up for matrix free algebraic combined method(ACM). Based on CTF phase 
% retrieval duality method for measured intensity data from one distance
%       
% In:   IR      intensity measurements
%       out     output from forward model PCT_forward
%       R       discrete CT projection matrix
%       sigma   propotionality constant: -delta/beta
%
% Out: 
%       A       system matrix function handle for ACM
%       b       data for ACM
%       AT      transposed system matrix function handle for ACM
%
% By Rasmus Dalgas Kongskov, 05/06/2014, DTU

% initialiez variables
[Nt,Nth] = size(IR);

% CTF model
P = 2*sigma*sin(out.w)-2*cos(out.w);

% mask to avoid near-zero crossing amplification
Z  = abs(P)>0;

% discrete measurements, dicrete fourier transform and masked
b = Z.*(1/sqrt(Nt)*fft(IR-1));

% model/system matrix set up
A  = @(x) Z.*P.*(1/sqrt(Nt)*fft(reshape(R*x(:),Nt,Nth)));
AT = @(x) reshape(R'*reshape(sqrt(Nt)*ifft((Z.*x).*P),Nt*Nth,1),out.N,out.N);

end