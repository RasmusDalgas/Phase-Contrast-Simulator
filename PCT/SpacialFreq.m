function [w,omega1,omega2] = SpacialFreq(ps,dim,lambda,R)
% In:   ps      pixelsize(m)
%       dim     object dimensions [n,m]
%       lambda  wavelength(m)
%       R       distance from object
% Out:  
%       w       pi*lambda*R*omega (omega = spacial frequencies)
%       omega   spacial frequencies (direction 1 or 2)
%
% By Rasmus Dalgas Kongskov, 08/03/2014, DTU

Fs = 1/ps; % sampling distance

% spacial frequencies
omega1 = Fs/dim(1):Fs/dim(1):Fs/2;
omega1 = [-fliplr(omega1),0,omega1(1:end-1)];
omega2 = Fs/dim(2):Fs/dim(2):Fs/2;
omega2 = [-fliplr(omega2),0,omega2(1:end-1)];

% values for phase retrieval
w1 = fftshift(pi*lambda*R*omega1.^2);
w2 = fftshift(pi*lambda*R*omega2.^2);
[W1,W2] = meshgrid(w2,w1);
w = W1+W2;

end