% Authors: R.Sacco, G.Guidoboni and A.G.Mauri
%
% This script implements the richardson stationary
% methods for linear solvers
%
function [x,err,niter] = richardson_stat(A,b,P,alpha,toll,nitmax)
% [x,err,niter] = richardson_stat(A,b,P,alpha,toll,nitmax)
err    = toll + 1 ;
niter  = 0 ;
n      = max(size(b)); 
x      = zeros(n,1);   
r      = b - A*x;
while ((err > toll) & (niter < nitmax))
    niter = niter + 1;
    z     = P \ r;
    x     = x + alpha*z;			
    r     = r - alpha*A*z;
    e     = norm(r)/norm(b);
    err   = [err; e];
end
