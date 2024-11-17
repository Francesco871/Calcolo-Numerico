%
close all
clear all
clc
%
a   = 0; 
b   = 1;
% 
n   = 3; 
h   = (b-a)/n; 
% interpolation nodes of the global polynomial interpolation
xn  = [a:h:b];
%
[l_i, x] = eval_global_l_i(a, b, n, xn);
%
setfonts;
%
N       = numel(x);
y_one   = ones(N,1);
y_zero  = zeros(N,1); 
yn_zero = zeros(n+1,1);
%
plot(x, l_i, x, y_one, 'k--', x, y_zero, 'k-', xn, yn_zero, 'om')
%
xlabel('x')
title('l_i(x)')
%
return
