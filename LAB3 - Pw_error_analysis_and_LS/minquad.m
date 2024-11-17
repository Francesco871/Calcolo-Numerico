% this script computes the least squares approximation of degree n of the m+1 pairs of data (x_i, y_i)
close all
clear all
clc
%
x      = [-pi:0.0001:pi]';
y      = sin(x)+cos(x.^2);
%
n      = 3;
%
[a, B] = least_squares(x, y, n);
%
c      = flip(a);
pstar  = polyval(c, x);
%
setfonts;
plot(x, y, x, pstar);
%  
return
