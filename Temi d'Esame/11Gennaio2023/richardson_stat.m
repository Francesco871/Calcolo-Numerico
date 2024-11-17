function [x,err,niter] = richardson_stat(A,b,P,alpha,toll,nitmax)
% Authors: R.Sacco
%
% This function implements the Richardson stationary iterative method for linear solvers
%
% [x,err,niter] = richardson_stat(A,b,P,alpha,toll,nitmax)

% initialization
err    = toll + 1 ;
niter  = 0 ;
n      = max(size(b)); 
x      = zeros(n,1);   
r      = b - A*x;
% iteration
while ((err > toll) & (niter < nitmax))
    niter = niter + 1;
    z     = P \ r;
    x     = x + alpha*z;			
    r     = r - alpha*A*z;
    e     = norm(r)/norm(b);
    err   = [err; e];
end
%
return
