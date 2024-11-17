% This script uses the function polyfit and polyval to evaluate the global 
% polynomial interpolation of degree n of the Runge function f(x) = 1/(1+x^2) 
% on the interval [a,b]=[-5,5]. The degree n is an integer >=1.
% Tchebicev non equally spaced nodes are used.
clear all
close all
clc
% define the interval 
xa    = -5;
xb    = +5;
% degree of the polynomial
n     = 10;
% Tchebicev nodes
theta = pi/n;
xi    = (xa+xb)/2 - (xb-xa)/2*cos([0:n]*theta);   
% define the function
fun   = @(x) 1./(1+x.^2);
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
legend('$x_i$ ', '$f(x)=1/(1+x^2)$', '$\Pi_nf(x)$','interpreter','latex')
%
grid
%
return
