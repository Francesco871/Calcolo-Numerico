% This script computes the integral of a function f=f(x) over the interval [a,b] and evaluates the error
clear all
close all
clc
%
a   = 0;
b   = +pi;
f   = @(x) sin(x) + cos(x);
Iex = -cos(b) + cos(a) + sin(b) - sin(a); 
%
% select quadrature formula:
%
% form = 'MP' -----> midpoint
% form = 'TR' -----> trapezoidal 
% form = 'CS' -----> Cavalieri-Simpson
% form = 'GL' -----> Gauss-Legendre
form = 'GL';
%
switch form
  %
  case 'MP' 
       % weights and nodes of the midpoint quadrature rule over [-1, +1]
       wn = 2;
       xn = 0;
  case 'TR' 
       % weights and nodes of the trapezoidal quadrature rule over [-1, +1]
       wn = [1;  1];
       xn = [-1; +1];
  case 'CS' 
       % weights and nodes of the Cavalieri-Simpson quadrature rule over [-1, +1]
       wn = 2*[1/6; 4/6; 1/6];
       xn = [-1; 0; +1];
  case 'GL' 
       % weights and nodes of the Gauss-Legendre quadrature rule over [-1, +1]
       wn = [1; 1];
       xn = [-1/sqrt(3); +1/sqrt(3)];
end
%
Mvec = [20:20:1000];
%
for ie = 1:numel(Mvec)
    Mh      = Mvec(ie);
    hi      = (b-a)/Mh;
    xv      = [a:hi:b];
    xb      = (xv(2:end) + xv(1:end-1))/2;
    H       = diff(xv);
    % evaluate the approximation of the integral
    Q       = quadrature(xb, H, f, wn, xn);
    % compute the absolute value of the error
    err(ie) = abs(Iex - Q);
    %
    h(ie)   = hi;
    %
    clear xv xb
end
%
setfonts
loglog(h, err, 'r-')
xlabel('$\log(h)$','interpreter','latex')
ylabel('$\log(\textrm{error})$','interpreter','latex')
%
return
