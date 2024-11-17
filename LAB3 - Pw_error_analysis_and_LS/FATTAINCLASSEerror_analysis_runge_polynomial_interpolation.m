% This script performs the error analysis of the interpolation error associated with 
% pw polynomial interpolation of degree r >=0. 
% The function f(x) = sin(omega*x) is studied on the interval [a,b]=[-pi,pi]. 
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
Mvec  = [1:1000];
% polynomial degree r
r     = 1;
% error
pw_interpolation_error = zeros(size(Mvec));
% loop over the (increasing) number of elements of the partition
for ie = 1:numel(Mvec)
    % number of elements of Tau_H(i)
    Mh    = Mvec(ie);
    % uniform mesh size
    h     = (xb-xa)/Mh; 
    H(ie) = h;
    % vertices of Tau_H(i)
    xv    = [xa:h:xb]';
    % distance between interpolation nodes over each element of Tau_H(i)
    if (r~=0)
       delta = h/r;
    end   
    % loop over each element of the partition Tau_H(i)
    err_old = 0;
    %
    for i = 1:Mh
        % endpoints of K_i
        xi    = xv(i);
        xip1  = xv(i+1);
        % pw interpolation nodes over K_i
        if (r~=0)
           xn = [xi:delta:xip1]';
        else
           xn = (xi+xip1)/2;
        end
        % nodal values of f over K_i
        yn    = fun(xn); 
        % evaluation of the r+1 local basis function set
        dx    = h/1000;
        % fine grid for local error evaluation
        xx    = [xi:dx:xip1]';
        Nxx   = numel(xx);
        %
        if (r~=0)
            phi = eval_local_FE_basis_functions(xi, xip1, dx, r);
        else
            phi = ones(Nxx,1);
        end
        % evaluation of the pw interpolation polynomial Pi_h^rf(x) over K_i
        pol   = zeros(Nxx, 1);
        for n = 1:numel(xn)
            pol = pol + yn(n)*phi(:,n);
        end    
        % evaluate f(x) over each K_i
        ff    = fun(xx); 
        % evaluate the interpolation error on K_i
        err   = norm(ff - pol, 'inf');
        % evaluate the maximum error between the previous and the present interval
        err_old = max(err_old, err);
    end
    % interpolation error on the partition Tau_H(i)
    pw_interpolation_error(ie) = err_old;
end
% plots
setfonts;
%
loglog(H, pw_interpolation_error, 'ko', H, pw_interpolation_error, 'r-')
xlabel('$h$','interpreter','latex')
legend('$\| f - \Pi_h^r f \|_{C^0([a,b])}$','interpreter','latex')
%
grid
%
return
