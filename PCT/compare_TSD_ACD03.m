close all; clear; clc;
addpath('../AIRtools','../Projector-Pack-0.2');

%% Compare two stage method (TSD) with algebraic combined model (ACD)

% phantom and settings
Nx  = 200; rand('state',5);
im = phantomgallery('grains',Nx,35);
im(im<1/5) = 3;
im(im<2/5) = 2;
im(im<2)   = 1;
im = padarray(im,[Nx/2,Nx/2]); N = size(im,1); % insert "frame"

% parameters and indexes
param = [inf,0.5,1e-6,40e3]; % [r1,r2,resolution,energy]
index{1} = [1.6407775512258e-7,8.4269742031797E-12]; % polycarbonate 40 keV
index{2} = [3.0141182422983e-7,2.6761112204583e-10]; % silicon 40 keV
index{3} = [3.3703495679786e-7,2.3169060965266e-10]; % aluminium 40 keV

% calculate Fresnell region
lambda = 299792458/(param(4)/4.135667e-15); (param(3))^2/lambda;

%% forward problem

% tomography settings
tomo.theta = 0:0.5:179.5;
tomo.p     = round(5+sqrt(2)*N)+mod(round(5+sqrt(2)*N),2);

% forward model PCT
[IR,out,mu,phi,R] = PCT_forward(im,param,index,tomo);

% CT data forward problem
I  = mu;

%% add noise

N0 = 1e5; % number of photons

IRn = -log(poissrnd(exp(-IR)*N0)/N0);

%% phase retrieval and set up
%  for two stage (TS) method and setup for algebraic combined (AC) method

sigma2 = -index{2}(1)/index{2}(2);
I_CTF = CTF_Dual_1R(IRn,out,sigma2);
[A_ACD,b,AT_ACD] = CTF_Dual_ACM_mf(IRn,out,R,sigma2);

%% reconstructions

% function handle expressions
A_TSD  = @(x) reshape(R*x(:),size(IR,1),size(IR,2));
AT_TSD = @(x) reshape(R'*x(:),N,N);

% TV-regularization
Ni  = 20000;
tol = 1e-6;
alpha_TSD = 0.12;
alpha_ACD = 112;
K_TSD = 5e3;
K_ACD = 260;

tic
[x_TSD,Pn2_TSD,Dn2_TSD,Pn1_TSD,Dn1_TSD,K_TSD] = ...
                  CP_TV_Reg(A_TSD,AT_TSD,I_CTF,alpha_TSD,[N,N],Ni,tol,K_TSD);
time1=toc; tic;
[x_ACD,Pn2_ACD,Dn2_ACD,Pn1_ACD,Dn1_ACD,K_ACD] = ...
                      CP_TV_Reg(A_ACD,AT_ACD,b,alpha_ACD,[N,N],Ni,tol,K_ACD);
time2=toc;

id = (N-Nx)/2+1:(N+Nx)/2; % without frame id

save(sprintf('results/compare_TSD_ACD03_%d',1),...
        'x_ACD','Pn2_ACD','Dn2_ACD','Pn1_ACD','Dn1_ACD','K_ACD', ...
        'x_TSD','Pn2_TSD','Dn2_TSD','Pn1_TSD','Dn1_TSD','K_TSD', ...
        'param','lambda','alpha_ACD','alpha_TSD', 'index','im','id',...
        'out','Nx','Ni','N0','N','I_CTF','I','IR','b','IRn','Ni','tol',...
        'time1','time2')

%% visualize
if false
fa = lambda/(2*pi); % scale factor to get indexes beta and delta
c = [fa*min(min(out.mu_real(id,id)))*0.9,fa*max(max(out.mu_real(id,id)))*1.1];

figure(1)
subplot(121); imshowBW(fa*x_TSD(id,id),c)
subplot(122); imshowBW(fa*x_ACD(id,id),c)

e_TSD = norm(x_TSD(id,id)-out.mu_real(id,id))/norm(out.mu_real(id,id))
e_ACD = norm(x_ACD(id,id)-out.mu_real(id,id))/norm(out.mu_real(id,id))

figure(2);
subplot(121); loglog((Pn2_TSD.^2)+(Dn2_TSD.^2));
subplot(122); loglog((Pn2_ACD.^2)+(Dn2_ACD.^2));
end