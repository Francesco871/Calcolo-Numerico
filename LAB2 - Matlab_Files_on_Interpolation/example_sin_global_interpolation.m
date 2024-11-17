% This script uses the function polyfit and polyval to evaluate the global 
% polynomial interpolation of degree n of the function f(x) = sin(omega*x) 
% on the interval [a,b]=[-pi,pi]. The degree n is an integer >=1.
clear all
close all
clc
% define the interval 
xa    = -pi;
xb    = +pi;
% degree of the polynomial
n     = 5;
% equally spaced nodes
dx    = (xb-xa)/n;
xi    = [xa:dx:xb];   
%xi = linspace(xa,xb,n+1);
% define the function
omega = 3;
fun   = @(x) sin(omega*x);
% nodal values of the function
yi    = fun(xi);
% global polynomial interpolation
pn    = polyfit(xi, yi, n);
% fine grid to plot the function and the polynomial 
xx    = [xa:(xb-xa)/1000:xb];
% evaluation of the function on the fine grid
f      = fun(xx);
% evaluation of the polynomial on the fine grid
Pi_n_f = polyval(pn, xx);
% interpolation error
format short e
err   = norm(f - Pi_n_f, 'inf')
% plots
setfonts;
%
plot(xi, yi, 'ko', xx, f, 'r-', xx, Pi_n_f, 'b--')
xlabel('$x$','interpreter','latex')
legend('$f(x_i)$ ', '$\sin(\omega \, x)$', '$\Pi_nf(x)$','interpreter','latex')
%
grid
%
return
