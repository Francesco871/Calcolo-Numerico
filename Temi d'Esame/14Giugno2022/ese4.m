clc
close all
clear all

format long e;

a = 4;
b = 6;
Int = [4 , 6];
alfa = 5;
g = @(x) sinh(x-alfa);
g_prime = @(x) cosh(x-alfa);

q = (g(b)-g(a))/(b-a);
tfun = @(x) x - q^(-1)*g(x);

itmax = 100;
toll = 1e-10;
x0 = b;

[alfa_c, niter_fix , err_fix] =fixedpoint(x0 ,tfun,toll ,itmax);
%errore effettivamente commesso:
e_alfa = alfa - alfa_c(end);

return