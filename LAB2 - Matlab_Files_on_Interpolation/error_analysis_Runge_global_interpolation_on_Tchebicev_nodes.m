% This script performs the error analysis of the global interpolation error as a function of n.
% The functions polyfit and polyval are used to evaluate the global 
% polynomial interpolation of degree n of the Runge function f(x) = 1/(1+x^2) 
% on the interval [a,b]=[-5,5]. The degree n is an integer >=1.
% Tchebicev non equally spaced nodes are used.
clear all
close all
clc
% define the interval 
xa    = -5;
xb    = +5;
% define the function
fun   = @(x) 1./(1+x.^2);
% fine grid to plot the function and the polynomial 
xx    = [xa:(xb-xa)/1000:xb];
% evaluation of the function on the fine grid
f     = fun(xx);
% vector containing the (increasing) polynomial degree n
Nvec  = [1:40];
% loop over the polynomial degree
for i = 1:numel(Nvec)
    n  = Nvec(i);
    % Tchebicev nodes
    theta = pi/n;
    xi    = (xa+xb)/2 - (xb-xa)/2*cos([0:n]*theta);   
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
