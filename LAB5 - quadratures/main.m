% This script computes the numerical approximation of the integral of [a,b] of a given function f=f(x)
close all
clear all
clc
% define the data
a   = 0;
b   = pi;
f   = @(x) exp(x);
% evaluate (if possible) the integral of f over [a,b]
Iex = exp(pi) - 1;
% define the partition
Mh = 10;
h  = (b-a)/Mh;
xv = [a:h:b];
xb = (xv(2:end) + xv(1:end-1))/2;
H  = diff(xv);
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
% evaluate the approximation of the integral
Q  = quadrature(xb, H, f, wn, xn);
% compute the absolute value of the error
err = abs(Iex - Q),
%
return
