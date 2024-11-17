% This script performs the error analysis of the interpolation error associated with 
% pw polynomial interpolation of degree r >=0. 
% The function f(x) = sin(omega*x) is studied on the interval [a,b]=[-pi,pi], with omega = 5.
clear all
close all
clc
% define the interval 
xa    = -pi;
xb    = +pi;
% define the function
omega = 5;
fun   = @(x) sin(omega*x); %+ cos(omega*x.^2);
% fine grid to evaluate the function, the pw interpolating polynomial and the error over the reference element [0, 1]
csi   = [0:1/10000:1];
% vector containing the (increasing) number of elements of the partition
Mvec  = [100:100:2000];
% polynomial degree r
r     = 1;
% error
MM    = numel(Mvec);
pw_interpolation_error = zeros(MM, 1);
% loop over the (increasing) number of elements of the partition
for ie = 1 : MM
    % number of elements of Tau_H(i)
    Mh    = Mvec(ie);
    % uniform mesh size
    h     = (xb-xa)/Mh; 
    H(ie) = h;
    % vertices of Tau_H(i)
    xv    = [xa:h:xb]';
    % loop over each element of the partition Tau_H(i)
    err_old = 0;
    %
    for i = 1 : Mh
        % endpoints of K_i
        xi    = xv(i);
        xip1  = xv(i+1);
        % evaluate f(x) over each K_i. To this purpose, we use the affine map x = xi + h*xi
        X     = xi + h*csi;
        X     = X';
        ff    = fun(X); 
        % pw interpolation nodes over K_i
        if (r~=0)
           % distance between interpolation nodes over each element of Tau_H(i)
           delta = (xip1 - xi)/r;
           xn    = [xi : delta : xip1]'; 
        else
           xn = (xi + xip1)/2;
        end
        %
        Nxn   = numel(xn);
        % nodal values of f over K_i
        yn    = fun(xn);
        % evaluation of the r+1 local basis function set
        NX    = numel(X);
        %
        if (r~=0)
            phi = eval_local_FE_basis_functions(xn, X, r);
        else
            phi = ones(NX,1);
        end
        % evaluation of the pw interpolation polynomial Pi_h^rf(x) over K_i
        pol   = zeros(NX, 1);
        for n = 1 : Nxn
            pol = pol + yn(n)*phi(:,n);
        end
        % evaluate the interpolation error on K_i
        err     = norm(ff - pol, 'inf');
        % evaluate the maximum error between the previous and the present interval
        err_old = max(err_old, err);
    end
    % interpolation error on the partition Tau_H(i)
    pw_interpolation_error(ie) = err_old;
    % clear variables for further use
    clear X xv xn
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
