% This script uses the function interp1 to evaluate the pw polynomial interpolation of degree 1 
% of the Runge function f(x) = 1/(1+x^2) on the interval [a,b]=[-5,5]. 
clear all
close all
clc
% define the interval 
xa    = -5;
xb    = +5;
% number of elements
Mh    = 10;
% vertices of the partition
dx    = (xb-xa)/Mh;
xi    = [xa:dx:xb];   
% define the function
fun   = @(x) 1./(1+x.^2);
% fine grid to plot the function and to evaluate the pw linear interpolating polynomial 
xx    = [xa:(xb-xa)/1000:xb];
% nodal values of the function
yi    = fun(xi);
% pw polynomial interpolation of degree 1
Pi_h_1_f = interp1(xi, yi, xx);
% evaluation of the function on the fine grid
f        = fun(xx);
% interpolation error
format short e
err      = norm(f - Pi_h_1_f, 'inf')
% plots
setfonts;
%
plot(xi, yi, 'ko', xx, f, 'r-', xx, Pi_h_1_f, 'b--')
xlabel('$x$','interpreter','latex')
legend('$x_i$ ', '$f(x)=1/(1+x^2)$', '$\Pi_h^1 f(x)$','interpreter','latex')
%
grid
%
return
