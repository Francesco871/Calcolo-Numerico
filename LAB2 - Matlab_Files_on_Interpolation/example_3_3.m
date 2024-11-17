% this script implements Example 3.3
close all
clear all
clc
%
fun   = @(x) sin(x);
%
a     = 0;
b     = 10;
%
nvect = [1:30];
%
err   = [];
xx    = [a:(b-a)/10000:b];
fx    = fun(xx);
%
for i=1:numel(nvect)
    %
    n    = nvect(i);
    xi   = [a:(b-a)/n:b];
    %
    X    = vander(xi);
    K(i) = cond(X);
    yi   = fun(xi);
    %
    pn   = polyfit(xi, yi, n);
    %
    pp   = polyval(pn, xx);
    %
    Enf  = fx - pp;
    %
    err  = [err; norm(Enf, 'inf')];		
end
%
setfonts;
%
figure;
semilogy(nvect, err, 'k-')
xlabel('n')
%
figure;
semilogy(nvect, K, 'ok--')
xlabel('n')
%
return
