clc
close all
clear all

T = 50;
MT = 1000;
t0 = 0;
Tfinal = t0+T;

f = @(t,y) -2.*t.*y.^(2);
dfdy = @(t,y) -4.*t.*y;

yex = @(t) 1./(1+t.^(2));
y0 = yex(0);

%1-2. Eulero in Avanti
thetaFE = 0;
[t,uFE] = theta_method(thetaFE, t0, Tfinal, y0, MT, f, dfdy);

%3-4. Eulero all'Indietro
thetaBE = 1;
[t,uBE] = theta_method(thetaBE, t0, Tfinal, y0, MT, f, dfdy);

%5-6. di Crank-Nicols
thetaCN = 1/2;
[t,uCN] = theta_method(thetaCN, t0, Tfinal, y0, MT, f, dfdy);

%7. yex
yy = yex(t);

setfonts;
figure
plot(t, yy, 'k-', t, uFE, 'b-', t, uBE, 'r-', t, uCN, 'm-');
legend('y(t)', 'u_{FE}', 'u_{BE}', 'u_{CN}')

%8. Errore FE
max_err_FE = norm(yy-uFE,'inf');

%9. Errore BE
max_err_BE = norm(yy-uBE,'inf');

%10. Errore CN
max_err_CN = norm(yy-uCN,'inf');

[max_err_FE , max_err_BE , max_err_CN]

return