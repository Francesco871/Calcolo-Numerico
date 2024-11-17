% This script computes the derivative of the function
% y(t) = sin(t) at t=t_0=pi/4 using formula (3.27e)
%
% x_h = (y(t_0+h) - y(t_0))/h        (*)
%
% and computes the error as a function of h, showing 
% first-order convergence of the method as h ----> 0.
close all
clear all
clc
%
setfonts;
h0  = 0.1;
h   = h0*[1; 1/2; 1/4; 1/8; 1/16; 1/32; 1/64; 1/128; 1/256];
%definition of the function
y         = @(t) sin(t);
%definition of the exact derivative of the function
yprime_ex = @(t) cos(t);
%
t0        = pi/4;
% compute the approximate derivative using formula (*)
yprime_h  = (y(t0+h) - y(t0))./(h);
%
format short e
% evaluate the error
err = yprime_ex(t0) - yprime_h

% plot the error as a function of h in log-log scale (cf. Figure 3.3)

figure;
loglog(h, abs(err))
xlabel('log(h)')
ylabel('log(err)')
%
return
