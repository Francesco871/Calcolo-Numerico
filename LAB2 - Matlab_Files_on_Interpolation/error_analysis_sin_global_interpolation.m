% This script performs the error analysis of the interpolation error in the global case as a function of n.
% The functions polyfit and polyval are used to evaluate the global 
% polynomial interpolation of degree n of the function f(x) = sin(omega*x) 
% on the interval [a,b]=[-pi,pi]. The degree n is an integer >=1.
clear all
close all
clc
% define the interval 
xa    = -pi;
xb    = +pi;
% define the function
omega = 3;
fun   = @(x) sin(omega*x);
% fine grid to plot the function and the polynomial 
xx    = [xa:(xb-xa)/1000:xb];
% evaluation of the function on the fine grid
f     = fun(xx);
% vector containing the (increasing) polynomial degree n
Nvec  = [1:20];
% loop over the polynomial degree
for i = 1:numel(Nvec)
    n  = Nvec(i);
    % equally spaced nodes
    dx = (xb-xa)/n;
    xi = [xa:dx:xb];   
    % nodal values of the function
    yi = fun(xi);
    % global polynomial interpolation
    pn = polyfit(xi, yi, n);
    % evaluation of the polynomial on the fine grid
    Pi_n_f = polyval(pn, xx);
    % interpolation error
    err(i) = norm(f - Pi_n_f, 'inf');
end
% plot of the error
setfonts;
%
semilogy(Nvec, err, 'ko', Nvec, err, 'b--')
xlabel('$n$','interpreter','latex')
title('$\| f - \Pi_n f \|_{C^0([a,b])}$','interpreter','latex')
%
grid
%
return
