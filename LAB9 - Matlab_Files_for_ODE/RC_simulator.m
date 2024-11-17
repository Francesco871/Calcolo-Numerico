% This script solves the Cauchy problem for a RC parallel circuit driven by an ideal current source.
% Two cases are considered: 
%
% 1: linear    resistor
% 2: voltage-controlled nonlinear resistor (pn junction diode)
%
clear all
close all
clc
% time interval []=s
t0     = 0;
T      = 100e-7;
tfin   = t0+T;
% input current []=A
ibar   = 100e-6;
is     = @(t) ibar.*(t > t0);
% capacitance []=F
C      = 1e-9;
% choice of the model for the resistor
imodel = 1;
%
if (imodel == 1)
   % linear resistor
   % resistance []=Ohm
   R      = 1e3;
   % conductance []=Siemens
   G      = 1/R;
   % characteristic of the resistor
   iR     = @(t,v) G*v + 0.*t + 0.*v;
else
   % nonlinear resistor: pn junction diode
   % reverse saturation current []=A
   i0     = ibar*1e-6;
   % Boltzmann's constant []=J/K
   kB     = 1.3806488e-23; 
   % elementary charge []=Coulomb
   q      = 1.6021892e-19;
   % temperature []=K
   temp   = 300;
   % thermal voltage
   Vth    = kB*temp/q;
   % characteristic of the resistor
   iR     = @(t,v) i0*(exp(v/Vth) - 1) + 0.*t;
end
% right-hand side of the Cauchy problem
RC_fun = @(t,v) 1/C*(is(t) - iR(t,v));
% initial condition: voltage across the capacitance []=V
v0     = 0;
%%%%%%%%%%%%%%%%%%
% MATLAB SOLVER
%%%%%%%%%%%%%%%%%%
tspan    = [t0, tfin];
[tM, vM] = ode15s(@(t,v) RC_fun(t, v), tspan, v0);
%%%%%%%%%%%%%%%%%%
% theta-method
%%%%%%%%%%%%%%%%%%
NT     = 300;
theta  = 0.5;
%
if (imodel == 1)
   dfdv = @(t, v) -G/C + 0.*t + 0.*v;
else
   dfdv = @(t,v) -i0/(C*Vth)*exp(v/Vth) + 0.*t;
end
%
[t, v] = theta_method(theta, t0, tfin, v0, NT, RC_fun, dfdv);
%
setfonts;
% set the axes
XMIN   = min([min(tM), min(t)]);
XMAX   = max([max(tM), max(t)]);
YMIN   = min([min(vM), min(v)]);
YMAX   = max([max(vM), max(v)]);
%
figure;
subplot(1,2,1)
plot(tM, vM, 'r-', tM, vM, 'ko')
axis([XMIN XMAX YMIN YMAX]);
xlabel('t [s]')
ylabel('v(t) [V]')
title('ode15s')
%
subplot(1,2,2)
plot(t, v, 'b-', t, v, 'ko')
axis([XMIN XMAX YMIN YMAX]);
xlabel('t [s]')
ylabel('v(t) [V]')
title('theta-metod')
%%%%%%%%%%%%%%%%%%
% ERROR ANALYSIS
%%%%%%%%%%%%%%%%%%
if (imodel == 1)
   % we can compute the exact solution
   tau    = R*C;
   vex    = @(t) R*ibar*(1 - exp(-(t-t0)/tau)) + v0*exp(-(t-t0)/tau);
   % errors
   eM     = vex(tM) - vM;
   etheta = vex(t)  - v;
   YMIN   = min([min(eM), min(etheta)]);
   YMAX   = max([max(eM), max(etheta)]);
   %
   figure
   subplot(1,2,1) 
   plot(tM, eM, 'r-', tM, eM, 'ko')
   axis([XMIN XMAX YMIN YMAX]);
   xlabel('t [s]')
   ylabel('error(t) [V]')
   title('ode15s')
   %
   subplot(1,2,2)
   plot(t, etheta, 'b-', t, etheta, 'ko')
   axis([XMIN XMAX YMIN YMAX]);
   xlabel('t [s]')
   ylabel('error(t)[V]')
   title('theta-metod')
   %
   [EM, index]     = max(abs(eM));
   time_M          = tM(index);
   [Etheta, index] = max(abs(etheta));
   time_theta      = t(index);
   EM, time_M 
   Etheta, time_theta
   %
end
%
return
