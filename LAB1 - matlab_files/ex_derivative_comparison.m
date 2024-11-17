% This script computes the derivative of the function
% y(t) = sin(t) at t=t_0=pi/4 using formula 
%
% yprimeh_forward =  (y(t_0+h) - y(t_0))/h        (*)
% yprimeh_backward = (y(t_0) - y(t_0-h))/h        (*)
% yprimeh_centered = (y(t_0+h) - y(t_0-h))/(2*h)        (*)
%
% and computes the error as a function of h, showing 
% the orders convergence of the methods as h ----> 0.
% N.B.: before running the script use the graphical command set path
% to add to the Matlab path the folder ./matlab_library


close all
clear all
%
setfonts;
h0  = 0.1;
h   = h0*[1; 1/2; 1/4; 1/8; 1/16; 1/32; 1/64; 1/128; 1/256];
%
y         = @(t) sin(t);
yprime_ex = @(t) cos(t);
t0        = pi/4;
% compute the approximate derivative using formula (*)
yprimeh_forward   = (y(t0 + h) - y(t0 ))./(h);
yprimeh_backward  = (y(t0 ) - y(t0-h ))./(h);
yprimeh_centered  = (y(t0 + h) - y(t0 - h))./(2*h);
%
format short e
% evaluate the error
err_forward  = abs(yprime_ex(t0) - yprimeh_forward);
err_backward = abs(yprime_ex(t0) - yprimeh_backward);
err_centered = abs(yprime_ex(t0) - yprimeh_centered);
%
%compute the slope of the three curves 
slope_forward =(log(err_forward(1:end-1))-log(err_forward(2:end)))./(log(h(1:end-1))-log(h(2:end)))
slope_backward=(log(err_backward(1:end-1))-log(err_backward(2:end)))./(log(h(1:end-1))-log(h(2:end)))
slope_centered=(log(err_centered(1:end-1))-log(err_centered(2:end)))./(log(h(1:end-1))-log(h(2:end)))
%
% plot the error as a function of h in log-log scale (cf. Figure 3.3)
setfonts;
figure;
loglog(h, err_forward,h, err_backward,h, err_centered )
xlabel('log(h)')
ylabel('log(err)')
legend('err forward','err backward','err centered') 
%
return