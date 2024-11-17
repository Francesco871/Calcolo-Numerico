% This script numerically solves the Cauchy problem
%
%      y^prime(t) = -K*y(t)*(M - y(t))    t>0
%      y(0)       = y_0
%
% using the theta method and the Matlab function ode15s
%
close all
clear all
clc
% data of the Cauchy problem
t_0     = 0;
T_final = 10;
y0      = 1;
M       = 20;
K       = 0.1;
A       = (M-y0)/y0;
%
yex     = @(t) M./(1 + A*exp(K*M*t));
%
f       = @(t,y) -K*y.*(M-y) + 0.*t;
dfdy    = @(t,y) -K*M + 2*K*y + 0.*t;
% evaluation of lambda
tt      = [t_0:(T_final-t_0)/10000:T_final];
%
Y       = yex(tt);
lambda  = norm(dfdy(tt,Y), 'inf')
% evaluation of dt_max for the Forward Euler method to be absolutely stable
dt_max  = 2/lambda
% determine the number of intervals for time discretization
% we set dt = K*dt_max and we compute NT by inverting the relation dt = (T_final - t_0)/NT
Kappa   = 0.1;
dt      = Kappa*dt_max;
NT      = round((T_final-t_0)/dt)
% choice of theta
theta   = 0;
[t, u]  = theta_method(theta, t_0, T_final, y0, NT, f, dfdy);
% use Matlab solver
[tM, uM] = ode15s(f, [t_0, T_final], y0);
%
setfonts;
%
figure
plot(tt, Y, 'm-', t, u, 'bo-', tM, uM, 'ko-')
xlabel('t')
legend('y(t)', '\theta method', 'ode15s')
%
return
