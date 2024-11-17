% Author: R.Sacco
%
% This script numerically computes the zero of g(x) = e^(-x) - x using:
%
% 1) graphical method
% 2) fzero
% 3) T(x) = exp(-x)
% 4) T(x) = Newton's method
% 
close all
clear all
clc

% g = g(x)
g     = @(x) exp(-x) - x;

% graphical preliminary study of g=g(x)
xplot = [-1:0.0001:1];

setfonts;
%
figure
plot(xplot, g(xplot), 'b-', xplot, zeros(size(xplot)), 'm-');
xlabel('x')

pause

% Matlab function fzero

% search interval
interval = [0 1];

% Define the options to use fzero (type help fzero for further indication)
options = optimset('Display','iter'); 

% Find the zero of g using fzero
[x_fzero, fval] = fzero(g, interval, options)

pause

% we set alpha = x_fzero for comparison with the other methods
alpha = x_fzero;

% fixed point iterations

% Parameters to monitor the convergence of the fixed-point iteration
toll  = 1e-12;
itmax = 10000;

% iteration function T = T(x)
T     = @(x) exp(-x);

x0    = 0;
% Find the zero of g using the fixed point iteration method
[x_fix, niter_fix, err_fix] = fixedpoint (x0, T, toll, itmax);

% evaluate the true error
alpha_fix = x_fix(end);
err_true  = alpha - alpha_fix;

% evaluate the estimated error
err_est   = err_fix(end);

niter_fix, alpha, alpha_fix, err_true, err_est

pause

% Newton's fixed point iteration

% iteration function T = T_N(x)
gpr   = @(x) -exp(-x) - 1;
T_N   = @(x) x - g(x)./gpr(x);

% Find the zero of g using Newton's method
[x_N, niter_N, err_N] = fixedpoint (x0, T_N, toll, itmax);

% evaluate the true error
alpha_N  = x_N(end);
err_true = alpha - alpha_N;

% evaluate the estimated error
err_est  = err_N(end);

niter_N, alpha, alpha_N, err_true, err_est

%
return
