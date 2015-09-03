function [IR,out,mu,phi,R] = PCT_forward(im,param,index,tomo)
% in:   im      2D phantom with values [0,1,2] 
%       param   vector of experimental parameters
%               [source to object (m), detector to object(m),...
%                pixelsize (m), energy (eV), sample-size(m)]
%       index   [refractive index decrement, arbsoroption index]
%       tomo    tomo.theta = angles(deg), tomo.p = number of projections
% Out: 
%       IR      Intensity image measured at detector
%       out     out.beta = absorption index, out.delta = decrement index, 
%               out.w = pi*lambda*R*omega (omega = spacial frequencies),
%               out.E0 = field just after object
%       mu      attenuation sinogram
%       phi     phase shift sinogram
%
% By Rasmus Dalgas Kongskov, 07/05/2014, DTU

% physical constants
h = 4.135667e-15;   % plancs constant (eV/s)
c = 299792458;      % speed of light (m/s)

% image and sinogram sizes
dim = size(im);
N   = dim(1);

if nargin<2
    % Default setup
    % simulated experiment settings
    r1 = inf;       % distance from source to object (inf = parallel beam)
    r2 = 1;         % distance from detector to object (m)
    ps = 1e-6;      % pixelsize (m)
    energy = 40e3;  % energy x-ray (eV)
else
    % simulated experiment settings
    r1 = param(1);      % distance from source to object
    r2 = param(2);      % distance from detector to object
    ps = param(3);      % pixelsize (m)
    energy = param(4);  % energy x-ray (eV)
end

if nargin<3
    % Default setup (non-realistic) indicies
    delta1 = 0;
    beta1  = 0;
    delta2 = 5e-9;
    beta2  = 5e-11;
    delta3 = 5e-9;
    beta3  = 5e-11;
else
    % indicies
    delta1 = index{1}(1);
    beta1  = index{1}(2);
    delta2 = index{2}(1);
    beta2  = index{2}(2);
    delta3 = index{3}(1);
    beta3  = index{3}(2);
end

if nargin<4
    % Defualt setup
    % parallel beam radon transform settings
    theta = 0:179; % angles
    Nt = round(5+sqrt(2)*N)+mod(round(5+sqrt(2)*N),2);
    Nth = 180;
else
    theta = tomo.theta;
    Nt = tomo.p;
    Nth = length(theta);
end

% wave-length of x-ray (m)
lambda = c/(energy/h);

% radon transform scaled to physical units
conf.x    = [-ps*N/2,ps*N/2]; conf.y     = [-ps*N/2,ps*N/2];
conf.dims = [N,N];            conf.r     = sqrt(2)*1.04*[-ps*N/2,ps*N/2];
conf.nr   = Nt;               conf.theta = theta;
R = parbeam(conf);

im2 = im;
im2(im==1) = 2*pi/lambda*(-delta1 + 1i*beta1);
im2(im==2) = 2*pi/lambda*(-delta2 + 1i*beta2);
im2(im==3) = 2*pi/lambda*(-delta3 + 1i*beta3);

sino = reshape(R*reshape(im2,N^2,1),Nt,Nth);

% phase and attenuation
phi  = real(sino);
mu   = imag(sino);

% field image after x-ray interaction with object
E0 = exp(-mu+1i*phi);

% spacial frequencies, only using 1D along line of x-ray
w = SpacialFreq(ps,[Nt,1],lambda,r2);
w = repmat(w(:,1),1,length(theta));

% Fresnel Propagator (fourier domain)
PR_f = exp(-1i*w);

% Intensity image at distance R
IR = abs(ifft( fft(E0).* PR_f ) ).^2;

% real attenuation of object (ground truth)
mu_real  =  2*pi/lambda*(beta1 *(im==1) + beta2 *(im==2) + beta3*(im==3) );
phi_real = -2*pi/lambda*(delta1*(im==1) + delta2*(im==2) + delta3*(im==3));

% output for further use
out.mu_real  = mu_real;  % ground truth
out.phi_real = phi_real; % ground truth
out.w = w;              % scaled frequencies

out.index = index;
out.theta = theta;
out.N = N;

end