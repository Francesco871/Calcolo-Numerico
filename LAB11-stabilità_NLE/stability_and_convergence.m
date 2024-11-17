% This script numerically solves the Cauchy problem
%
%      y^prime(t) = -lambda*y(t)     t > t_0
%      y(t_0)     = y_0
%
% using the theta method 
%
close all
clear all
clc
% data of the Cauchy problem
t_0     = 0;
T_final = 50;
%
lambda  = 1;
%
y0      = 1;
yex     = @(t) y0*exp(-lambda*t);
%
f       = @(t,y) -lambda*y + 0.*t;
dfdy    = @(t,y) -lambda + 0.*y + 0.*t;
% evaluation of lambda
tt      = [t_0:(T_final-t_0)/10000:T_final];
%
Y       = yex(tt);
% evaluation of dt_max for the Forward Euler method to be absolutely stable
dt_max  = 2/lambda;
%
K       = [0.1, 0.05, 0.025, 0.0125, 0.0125/2, 0.0125/4, 0.0125/8, 0.0125/16];
errore  = [];
deltat  = [];
% choice of theta
theta   = 0;
%
% convergence analysis as a function of h=dt
for i = 1:numel(K)
    % determine the number of intervals for time discretization
    % we set dt = K*dt_max and we compute NT by inverting the relation dt = (T_final - t_0)/NT
    Kappa   = K(i);
    dt      = Kappa*dt_max;
    deltat  = [deltat, dt];
    %
    NT      = round((T_final-t_0)/dt); 
    [t, u]  = theta_method(theta, t_0, T_final, y0, NT, f, dfdy);
    %
    errore  = [errore, norm(yex(t)-u, 'inf')];
%
end
setfonts;
%
loglog(deltat, errore, '*r-')
xlabel('delta t')
ylabel('errore')
% calcolo di p
p = log(errore(2:end)./errore(1:end-1))./log(deltat(2:end)./deltat(1:end-1))
%
return
