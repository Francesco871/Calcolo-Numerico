% Authors: R.Sacco
%
% This script performs an experimental study of the convergence of an error sequence.
% The sequence is generated by the numerical approximation of the quantity
%
%     x = y'(t_0) 
%
% where:
%
%     y   = y(t) = exp(l*t)
%     t_0 = 1
%    
% using centered, forward and backward finite differences

clear all
close all
setfonts
% 
l  = 0.1;
t0 = 1;
%
y  = @(t) exp(l*t);
yp = @(t) l*exp(l*t);
%
x  = yp(t0);
% definition of h
n  = 8;
h0 = 0.1;
h  = h0./2.^[0:n]; 
% definition of t0+h and t0-h for the evaluation of the finite difference approximation of y'(t_0)
t0_plu_h = t0 + h;
t0_min_h = t0 - h;
%
xh_CD    = (y(t0_plu_h) - y(t0_min_h))./(2.*h);
xh_FD    = (y(t0_plu_h) - y(t0))./h;
xh_BD    = (y(t0) - y(t0_min_h))./h;
% error evaluation as a function of h
err_CD   = x - xh_CD;
err_FD   = x - xh_FD;
err_BD   = x - xh_BD;

% Compute the asymptotic order of convergence and error constant for the three methods: 
[p_FD, C_FD] = order_estimate(h, err_FD);
[p_BD, C_BD] = order_estimate(h, err_BD);
[p_CD, C_CD] = order_estimate(h, err_CD);

% Plots

% log-log plot of the errors
figure
loglog(h, abs(err_BD), 'ko-', h, abs(err_FD), 'bo-', h, abs(err_CD), 'go-' )
xlabel('$\ln(h)$',  'Interpreter', 'Latex')
ylabel('$\ln(err)$','Interpreter', 'Latex')
legend('BD method', 'FD method', 'CD method')

% plot of the asymptotic order of convergence p
figure
plot(h(2:end), p_BD, 'ko-', h(2:end), p_FD, 'bo-', h(2:end), p_CD, 'go-' )
xlabel('$h$',  'Interpreter', 'Latex')
ylabel('$p$','Interpreter', 'Latex')
legend('BD method', 'FD method', 'CD method')

% plot of the asymptotic error constant
figure
plot(h(2:end), C_BD, 'ko-', h(2:end), C_FD, 'bo-', h(2:end), C_CD, 'go-' )
xlabel('$h$',  'Interpreter', 'Latex')
ylabel('$C$','Interpreter', 'Latex')
legend('BD method', 'FD method', 'CD method')

return
