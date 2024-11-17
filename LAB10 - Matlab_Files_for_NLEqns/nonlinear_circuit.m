% This script computes the working point of the circuit constituted by a current source
% driving a parallel between a resistor and a p-n diode
clear all
close all
clc

% physical constants
T   = 300;                             % ref. temperature value ( [] = K) 
kB  = 1.3806488e-23;                   % Boltzmann constant ([] = V A s K^{-1} = m^2 kg s^{-2} K^{-1})
q   = 1.6021892e-19;                   % elementary charge ([] = C)
Vth = kB*T/q ;                         % thermal voltage ([] = V)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test case # 1
% input parameters of the circuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diode saturation current ([]=A)
is  = 1e-3;
% current source ([]=A)
i0  = 2*is;
% resistance ([]=Ohm)
R   = 1e2;

a   = R*(is + i0);
b   = R*i0;

% element I-V relationships
g1  = @(v) v;
g2  = @(v) a - b*exp(v/Vth);
g   = @(v) g1(v) - g2(v);

% voltage interval ([]=V)
v1  = -0.1;
v2  = +0.1;

% 
vv  = [v1:(v2-v1)/10000:v2];

% set screen parameters for the plots
setfonts;

% graphical study of the circuit
figure
plot(vv, g1(vv), 'b-', vv, g2(vv), 'm-') 
xlabel('v [V]')
ylabel('[V]')
legend('g1(v)', 'g2(v)')
grid

figure
plot(vv, g(vv), 'k-', vv, zeros(size(vv)), 'r-') 
xlabel('v [V]')
ylabel('[V]')
legend('g(v)')
grid

% numerical evaluation of the working point using fzero
format long e
vstar = fzero(g, [v1 v2])
%

% numerical evaluation of the working point using fixed-point iterations
T1 = @(v) a - b*exp(v/Vth);
T2 = @(v) Vth*log(1 + is/i0 - v/(R*i0));

% we study the applicability of T1 using Ostrowski lemma
T1_prime = @(v) -b/Vth.*exp(v/Vth);
abs(T1_prime(vstar))

% we study the applicability of T2 using Ostrowski lemma
T2_prime = @(v) -Vth./(1 + is/i0 - v/(R*i0))*1/(R*i0);
abs(T2_prime(vstar))

% we use T2 for the numerical approximation of vstar
v0    = v2;
itmax = 1000; 
toll  = 1e-15;

[vfix, niter, err] = fixedpoint (v0, T2, toll, itmax);

% computed approximation
vstar_fix          = vfix(end)
% true error
true_error         = vstar - vstar_fix
% estimated error
est_error          = err(end)
niter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test case # 2
% input parameters of the circuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diode saturation current ([]=A)
is  = 1e-6;
% current source ([]=A)
i0  = 1e-12;
% resistance ([]=Ohm)
R   = 1e5;

a   = R*(is + i0);
b   = R*i0;

% element I-V relationships
g1  = @(v) v;
g2  = @(v) a - b*exp(v/Vth);
g   = @(v) g1(v) - g2(v);

% voltage interval ([]=V)
v1  = -0.01;
v2  = +0.4;

% 
vv  = [v1:(v2-v1)/10000:v2];

% set screen parameters for the plots
setfonts;

% graphical study of the circuit
figure
plot(vv, g1(vv), 'b-', vv, g2(vv), 'm-') 
xlabel('v [V]')
ylabel('[V]')
legend('g1(v)', 'g2(v)')
grid

figure
plot(vv, g(vv), 'k-', vv, zeros(size(vv)), 'r-') 
xlabel('v [V]')
ylabel('[V]')
legend('g(v)')
grid

% numerical evaluation of the working point using fzero
format long e
vstar = fzero(g, [v1 v2])
%

% numerical evaluation of the working point using fixed-point iterations
T1 = @(v) a - b*exp(v/Vth);
T2 = @(v) Vth*log(1 + is/i0 - v/(R*i0));

% we study the applicability of T1 using Ostrowski lemma
T1_prime = @(v) -b/Vth.*exp(v/Vth);
abs(T1_prime(vstar))

% we study the applicability of T2 using Ostrowski lemma
T2_prime = @(v) -Vth./(1 + is/i0 - v/(R*i0))*1/(R*i0);
abs(T2_prime(vstar))

% we use T1 for the numerical approximation of vstar
v0    = v2;
itmax = 1000; 
toll  = 1e-15;

[vfix, niter, err] = fixedpoint (v0, T1, toll, itmax);

% computed approximation
vstar_fix          = vfix(end)
% true error
true_error         = vstar - vstar_fix
% estimated error
est_error          = err(end)
niter

%
return
