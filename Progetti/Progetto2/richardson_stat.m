function [x,err,k] = richardson_stat(A,b,x0,P,alpha,toll,nitmax,stop_test)
%
% Authors: Riccardo Sacco
%
% This script implements the stationary Richardson method
%
% The input parameter stop_test is a flag to define the termination test of the algorithm:
%
%     stop_test = 1 ------> termination criterion based on the control of the increment
%     stop_test = 2 ------> termination criterion based on the control of the residual
%
% [x,err,k] = richardson_stat(A,b,x0,P,alpha,toll,nitmax,stop_test)
%
n      = max(size(b)); 
% initial guess
x      = x0;   
% initial residual
r      = b - A*x;
% initialization for the termination test
r0     = r;
xold   = x;
%
k      = 0 ;
e      = toll + 1;
err    = [];

%
while ((e > toll) & (k < nitmax))
    % step from k to k+1
    z       = P \ r;
    w       = A*z;
    x       = x + alpha*z;			
    r       = r - alpha*w;
    %
    if (stop_test == 1)
       % the convergence test is based on the control of the norm of the increment
       e    = norm(x - xold);
    else
       % the convergence test is based on the control of the norm of the residual
       e    = norm(inv(P)*r)/norm(inv(P)*r0);
    end
    %
    err   = [err; e];
    k     = k + 1;
    xold  = x;
    %
end
%
return
