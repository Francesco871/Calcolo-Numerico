% This script computes the derivative of the function
% y(t) = sin(t) at t=t_0=pi/4 using formula 
%
% x_h = (y(t_0+h) - y(t_0-h))/(2*h)        (*)
%
% and computes the error as a function of h, showing 
% second-order convergence of the method as h ----> 0.
% N.B.: before running the script use the graphical command set path
% to add to the Matlab path the folder ./matlab_library
close all
clear all
%
h0  = 0.1;
h   = h0*[1; 1/2; 1/4; 1/8; 1/16; 1/32; 1/64; 1/128; 1/256];
%
y         = @(t) sin(t);
yprime_ex = @(t) cos(t);
t0        = pi/4;
% compute the approximate derivative using formula (*)
yprime_h  = (y(t0 + h) - y(t0 - h))./(2*h);
%
format short e
% evaluate the error
err = yprime_ex(t0) - yprime_h
% plot the error as a function of h in log-log scale (cf. Figure 3.3)
setfonts;
figure;
loglog(h, abs(err))
xlabel('log(h)')
ylabel('log(err)')
%
return
