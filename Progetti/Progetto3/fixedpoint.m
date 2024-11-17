function [x, niter, err] = fixedpoint (x0, T, toll, itmax)
% Author: R. Sacco
%
% This function implements the fixed point iteration method to numerically compute the zero of g=g(x)
% 
% [x, niter, err] = fixedpoint (x0, T, toll, itmax)
%
x    = x0;
xold = x;
%
k    = 0;
e    = toll + 1;
err  = [];
%
while ((k < itmax) & (abs(e) > toll))
      xnew = T(xold);
      e    = xnew - xold;
      x    = [x; xnew];
      err  = [err; e];
      xold = xnew;
      k    = k+1;
end
%
niter = k;
%
return
