% This script verifies the necessary and sufficient condition to be satisfied by the acceleration parameter 
% alpha for the Stationary Richardson method to converge
clear all
close all
clc
% definition of the eigenvalues of P^{-1}*A
n  = 7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The case of complex eigenvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Re = [-7, -5, -3, -1, 2, 5, 6];
Im = 5*sin([1:n]);
mu = complex(Re, Im);
%
setfonts;
% plot of the eigenvalues in the complex plane
figure
plot(mu, 'o')
xlabel('Re(z)')
ylabel('Im(z)')
title('\{\mu_k\}')
axis('square')
%
index = [1:n];
f     = 2*real(mu)./abs(mu).^2;
uno   = ones(size(index)); 
% several choices of \alpha
alpha = 10;
ff    = f/alpha;
% verification of the theorem
ff > 1
%
figure
plot(index, ff, 'b--', index, ff, 'ro', index, uno, 'k-')
xlabel('eigenvalue index k')
ylabel('2*Re(\mu_k)/(\alpha|\mu_k|^2)')
%
alpha = -10;
ff    = f/alpha;
% verification of the theorem
ff > 1
%
figure
plot(index, ff, 'b--', index, ff, 'ro', index, uno, 'k-')
xlabel('eigenvalue index k')
ylabel('2*Re(\mu_k)/(\alpha|\mu_k|^2)')
%
alpha = 0.1;
ff    = f/alpha;
% verification of the theorem
ff > 1
%
figure
plot(index, ff, 'b--', index, ff, 'ro', index, uno, 'k-')
xlabel('eigenvalue index k')
ylabel('2*Re(\mu_k)/(\alpha|\mu_k|^2)')
%
alpha = -0.1;
ff    = f/alpha;
% verification of the theorem
ff > 1
%
figure
plot(index, ff, 'b--', index, ff, 'ro', index, uno, 'k-')
xlabel('eigenvalue index k')
ylabel('2*Re(\mu_k)/(\alpha|\mu_k|^2)')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The case of real eigenvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = Re;
f  = 2./mu;
%
alpha = 10;
ff    = f/alpha;
% verification of the theorem
ff > 1
%
figure
plot(index, ff, 'b--', index, ff, 'ro', index, uno, 'k-')
xlabel('eigenvalue index k')
ylabel('2/(\alpha \mu_k)')
%
alpha = -10;
ff    = f/alpha;
% verification of the theorem
ff > 1
%
figure
plot(index, ff, 'b--', index, ff, 'ro', index, uno, 'k-')
xlabel('eigenvalue index k')
ylabel('2/(\alpha \mu_k)')
%
alpha = 0.1;
ff    = f/alpha;
% verification of the theorem
ff > 1
%
figure
plot(index, ff, 'b--', index, ff, 'ro', index, uno, 'k-')
xlabel('eigenvalue index k')
ylabel('2/(\alpha \mu_k)')
%
alpha = -0.1;
ff    = f/alpha;
% verification of the theorem
ff > 1
%
figure
plot(index, ff, 'b--', index, ff, 'ro', index, uno, 'k-')
xlabel('eigenvalue index k')
ylabel('2/(\alpha \mu_k)')
%
return
