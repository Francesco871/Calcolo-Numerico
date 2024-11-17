clc
close all
clear all

setfonts;

%1)ANALISI CIRCUITALE IN REGIME STATICO: METODO GRAFICO
E = 1.5; %[V]
R1 = 100; %[ohm]
R2 = 100; %[ohm]
KB = 1.3806488e-23; %[J*K^-1]
q =  1.6021892e-19; %[C]
T = 300; %[K]
i0 = 1e-3; %[A]

%1.Definisco A, B, M e Vth
A = (E/R1)+i0;
B = (1/R1)+(1/R2);
M = i0;
Vth = KB*T/q;

%2.Plotto g1(v) e g2(v) e trovo vstar_grafico con il metodo grafico
v1 = -0.1; %[V]
v2 = +0.1; %[V]
J = [v1 , v2];
vv = [v1:(v2-v1)/10000:v2];

g1 = @(v) A -B.*v;
g2 = @(v) M*exp(v./Vth);

figure
plot(vv, g1(vv), 'b', vv, g2(vv), 'k')
xlabel('v[V]')
ylabel('g')
legend('g1' , 'g2')
grid
%vstar_grafico = 0.06933284

%2)ANALISI CIRCUITALE IN REGIME STATICO: METODO NUMERICO
format long e;

%1.Uso fzero per calcolare vstar_matlab
g = @(v) g1(v)-g2(v);
[vstar_matlab, gval] = fzero(g, J);

%2.Definisco T1'(v)
T1_prime = @(v) -R1*i0/(Vth*2)*exp(v./Vth);

%3.Calcolo y1
y1 = abs(T1_prime(vstar_matlab));
%Il metodo T1 non è convergente perché y1>1 (metodo di Ostrowski)

%4.Definisco T2'(v)
T2_prime = @(v) -2*Vth./(E+R1*i0-2.*v);

%5.Calcolo y2
y2 = abs(T2_prime(vstar_matlab));
%Il metodo T2 è convergente perché y2>1

%6.Eseguo la function fixedpoint per calcolare vstar_fix, calcolo poi
%  l'est_err e true_err
v0 = -0.1;
itmax = 1000;
toll = 1e-10;
T2 = @(v) Vth*log((E+R1*i0-2.*v)./(R1*i0));

[vfix, niter, err] = fixedpoint (v0, T2, toll, itmax);

vstar_fix = vfix(end);
est_err = err(end);
true_err = vstar_matlab - vstar_fix;

%3)ANALISI CIRCUITALE IN REGIME DINAMICO
C = 1e-3; %[F]
t0 = 0; %[s]
Tfinal = 0.6; %[s]
It = [t0 , Tfinal];

imax = 1e-2; %[A]
t1 = 0.1; %[s]
t2 = 0.3; %[s]
t3 = 0.5; %[s]

is = @(t) 0 -imax/3 * (t > t1 & t <= t2) +imax * (t > t2 & t <= t3);

%1.Definisco f(t,v)
f = @(t,v) is(t)/C - v/(R1*C) - (i0/C)*(exp(v/Vth)-1);

%2.Uso il comando theta_method per calcolare l'approssimazione numerica del
%  PC, plotto la soluzione e calcolo deltaV
theta = 0.5;
v0 = vstar_matlab;
NT = 2000;

dfdv = @(t,v) -1/(R1*C) - (i0/(C*Vth))*exp(v/Vth);

[t, V] = theta_method(theta, t0, Tfinal, v0, NT, f, dfdv);

figure
plot(t, V, 'm', t, 0*t, 'b')
xlabel('t[s]')
ylabel('v[V]')

deltaV = max(V) - min(V);

return