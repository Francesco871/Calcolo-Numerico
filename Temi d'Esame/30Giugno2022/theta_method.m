function [t,u] = theta_method(theta, t0, Tfinal, y0, NT, f, dfdy)
% Authors: R.Sacco
%
% [t,u] = theta_method(theta, t0, Tfinal, y0, NT, f, dfdy)
%
% This script implements the Theta-method for the general Cauchy problem:
%
% y'(t) = f(t,y(t)) ;
% y(t0) = y0 ;
%
% theta is in [0,1] (theta=0 FE, theta=1 BE, theta=0.5 CN)
% t0 is the initial integration time, Tfinal=t0+T is the final integration time
% f    = @(t,y) ..... handle function of the rhs 
% dfdy = @(t,y) ..... handle function of the derivative of f
%                     with respect to the variable y
%

% Tolerance for the Newton fixed point iteration
tol   = 1e-10;
% Max number of Newton iterations
maxit = 100;
% initial value of Cauchy problem (array needed to store results)
u     = [y0];
% time discretization
h     = (Tfinal-t0)/NT;
% initial time of the integration (array needed to store time)
t     = [t0];

% Time Integration using theta method with Netwon algos
for n = 0:NT-1

% Newton initialization of the new time level
err = 1;
it  = 0;

% Initialization of the Newton method with the solution of the previuous
% time step integration
un  = u(n+1);
x   = un;

% Calculation of the present time and the following time 
tn  = t0+n*h;
tnp = tn+h;

% Updating time array
t   = [t; tnp];

% Newton's iteration
while ((err > tol) & (it < maxit))
   
   fn  = f(tn,un); % evalution of the rhs at present time step
   fnp = f(tnp,x); % evalution of the rhs at the following time step
   
   % Calculation of the funtion and its derivative at which zero must be found
   % that is the SOLVE STEP
   F   = x - un - h*(1-theta)*fn - h*theta*fnp;
   Fpr = 1 - h*theta*dfdy(tnp,x);
   dx  = -F/Fpr;
   % Calculation of the error at the present Newton's iteration step
   err = abs(dx);
   % Updating the solution
   % That is the UPDATE STEP
   x   = x+dx;
    
   it = it +1 ;
   
end

% Updating the vector of the time integrated solution
u = [u; x];
 
end

return
