% This script numerically solves the Cauchy problem
%
%      y^prime(t) = f(t, y(t))      t > t_0
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
yex     = @(t) 1./(1+t.^2);
y0      = yex(0);
%
f       = @(t,y) -2*t.*y.^2;
dfdy    = @(t,y) -4*t.*y;
% evaluation of lambda
tt      = [t_0:(T_final-t_0)/10000:T_final];
%
Y       = yex(tt);
lambda  = norm(dfdy(tt,Y), 'inf');
% evaluation of dt_max for the Forward Euler method to be absolutely stable
dt_max  = 2/lambda;
% determine the number of intervals for time discretization
% we set dt = K*dt_max and we compute NT by inverting the relation dt = (T_final - t_0)/NT
Kappa   = 1;
dt      = Kappa*dt_max;
%
NT      = round((T_final-t_0)/dt);
% choice of theta
theta   = 0;
[t, u]  = theta_method(theta, t_0, T_final, y0, NT, f, dfdy);
%
setfonts;
%
plot(tt, Y, 'm-', t, u, 'b-')
xlabel('t')
legend('y(t)', 'u_h')
%
return
