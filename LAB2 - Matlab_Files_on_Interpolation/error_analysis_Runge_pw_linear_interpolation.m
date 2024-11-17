% This script performs the error analysis of the interpolation error in the pw linear case.
% The function interp1 is used to evaluate the pw polynomial interpolation of degree 1 
% of the Runge function f(x) = 1/(1+x^2) on the interval [a,b]=[-5,5]. 
clear all
close all
clc
% define the interval 
xa    = -5;
xb    = +5;
% define the function
fun   = @(x) 1./(1+x.^2);
% fine grid to to evaluate the function and the pw linear interpolating polynomial 
xx    = [xa:(xb-xa)/10000:xb];
% vector containing the (increasing) number of elements of the partition
Mvec  = 99*[1:20];
% loop over the number of elements of the partition
for i = 1:numel(Mvec)
    % number of elements
    Mh    = Mvec(i);
    h     = (xb-xa)/Mh;
    H(i)  = h;
    % vertices of the partition
    xi    = [xa:h:xb];   
    % nodal values of the function
    yi    = fun(xi);
    % pw polynomial interpolation of degree 1
    Pi_h_1_f = interp1(xi, yi, xx);
    % evaluation of the function on the fine grid
    f        = fun(xx);
    % interpolation error
    err(i)   = norm(f - Pi_h_1_f, 'inf');
end
% plots
setfonts;
%
loglog(H, err, 'ko', H, err, 'r-')
xlabel('$h$','interpreter','latex')
legend('$\| f - \Pi_h^1 f \_{C^0([a,b])}$','interpreter','latex')
%
grid
%
return
