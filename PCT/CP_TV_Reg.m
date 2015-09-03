function [u,Pn2,Dn2,Pn1,Dn1,K] = CP_TV_Reg(A,AT,b,lambda,dim,Ni,tol,K)
% Matrix-free version of Chambolle-Pock algorithm for solving the
% TV-regularization problem:
%       
% min_x ||A(u) - b||_2^2 + lambda*TV(x)
%
% in:   A       function handle which takes x in dimensions NxN and
%               returns the sinogram of size Nt x Nth
%       AT      function handle which return the transpose of A
%       b       data term of size Nt x Nth
%       lambda  regularization parameter (positive)
%       dim     dimesions of solution u
%       Ni      maximum number of iterations
%       tol     tolerance for stopping criteria
%
% Out: 
%       ub      final solution
%       en      change of solution for each iteration
%
% By Rasmus Dalgas Kongskov, 23/10/2014, DTU

% finite diffrence derivative approximation matrix
D = sparse(-eye(dim(1))+diag(ones(dim(1)-1,1),1));
D = [kron(eye(dim(1)),D);kron(D,eye(dim(1)))];

% power method for calculating L=||(A,D)||_2
xL = ones(dim(1),dim(2));
for n = 1:50
    xL = AT(A(xL))+reshape(D'*(D*xL(:)),dim(1),dim(2));
    xL = xL/norm(xL);
end

% power norms of operators
nA = norm(A(xL));
nD = norm(reshape(D*xL(:),dim(1),dim(2)*2));

% scaling to make norms of A and D of similar magnitude
if true
if nA<nD
    A  = @(x) K*A(x);
    AT = @(x) K*AT(x);
    b  = K*b;
    nA = norm(A(xL));
else
    D = K*D;
    nD = norm(reshape(D*xL(:),dim(1),dim(2)*2));
end
end

L  = sqrt(nA^2+nD^2);

% parameters for method
tau   = 1/L;
sigma = 1/L;
theta = 1;

% initialize variables for method
u    = zeros(dim(1),dim(2));
p    = b*0;
q    = [u(:)*0;u(:)*0];
ub   = u;
Aubn = A(ub);
ATp  = AT(p);

% adaptivity varibles
Pn1 = zeros(1,Ni); Pn2 = zeros(1,Ni);
Dn1 = zeros(1,Ni); Dn2 = zeros(1,Ni);
a = 0.5;
delta = 1.5;
eta = 0.95;
s = 200;

% iterative CP loop
for n = 1:Ni
    
    % use method
    pn   = (p+sigma*(Aubn-b))/(1+sigma);
    ATpn = AT(pn);
    
    prq  = q+sigma*D*ub(:);
    q    = lambda*prq./max(lambda,abs(prq));
    
    un   = u  - tau*(ATpn+reshape(D'*q,dim(1),dim(2)));
    ubn  = un + theta*(un-u);
    Aub  = Aubn;
    Aubn = A(ubn);
    
    % primal/dual residuals
    Pn = 1/tau*(ub-ubn)-(ATp-ATpn);
    Dn = 1/sigma*(p-pn)-(Aub-Aubn);
    Pn1(n) = norm(Pn,1);
    Pn2(n) = norm(Pn,2);
    Dn1(n) = norm(Dn,1);
    Dn2(n) = norm(Dn,2);
    
    % adapt parameters
    if Pn1(n) > s*delta*Dn1(n)
        tau   = tau/(1-a);
        sigma = sigma*(1-a);
        a = a*eta;
    end
    
    if Pn1(n) < s/delta*Dn1(n)
        tau   = tau*(1-a);
        sigma = sigma/(1-a);
        a = a*eta;
    end
    
    if (Pn2(n)^2 + Dn2(n)^2)/(Pn2(1)^2 + Dn2(1)^2) < tol
        break
    end
     
    % step to next n
    u   = un;
    ub  = ubn;
    p   = pn;
    ATp = ATpn;
    
end

end