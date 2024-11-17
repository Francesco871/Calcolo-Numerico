clc
close all
clear all

format long e;

alfa = 3;

itmax = 100;
toll = 1e-12;
x0 = 1;

%1.
tfun1 = @(x) (2/3).*x -(1/3)*alfa;

[x_fix1, niter1 , err1] =fixedpoint(x0 ,tfun1,toll ,itmax);

alfa_c1 = x_fix1(end);

e_alfa1 = alfa - alfa_c1;

%2.
tfun2 = @(x) alfa + x.*0;

[x_fix2, niter2 , err2] =fixedpoint(x0 ,tfun2,toll ,itmax);

alfa_c2 = x_fix2(end);

e_alfa2 = alfa - alfa_c2;

return