% this script implements Example 3.1.
% The same example is solved by using the least-squares approximation of degree n < m
clear all
close all
clc
%
a  = 0;
b  = 1;
% the number of data is equal to m+1
m  = 10;
xi = [a:(b-a)/m:b]';
yi = [0.1576, 0.9706, 0.9572, 0.4854, 0.8003, 0.1419, 0.4218, 0.9157, 0.7922, 0.9595, 0.6557]';
% global interpolation: n=m
format short e
n  = m;
pn = polyfit(xi,yi,n);
cn = fliplr(pn);
%  
for i=1:(n+1)
    X(:, i) = xi.^(i-1);
end
%
c  = X \ yi;
%
[cn', c]
%
xx = [a:0.001:b];
pp = polyval(pn,xx);
% least-squares approximation of degree n < m
n      = 5;
pnstar = polyfit(xi,yi,n);
ppstar = polyval(pnstar,xx);
%
setfonts;
%
figure
%
subplot(1,2,1)
plot(xi, yi, 'ok', xx, pp, 'm-')
xlabel('$x$', 'Interpreter','latex')
legend('$y_i$', '$\Pi_n f$','Interpreter','latex')
%
subplot(1,2,2)
plot(xi, yi, 'ok', xx, ppstar, 'm-')
xlabel('$x$', 'Interpreter','latex')
legend('$y_i$', '$p_n^\ast$','Interpreter','latex')
%
return
