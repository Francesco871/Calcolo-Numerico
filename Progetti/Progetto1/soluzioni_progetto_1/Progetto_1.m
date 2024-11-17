% Questo script contiene la soluzione del progetto #1
clear all
close all
format short e
%
setfonts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Analisi, implementazione e visualizzazione
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IS = 1e-6;     % A
R  = 2e3;      % Ohm
C  = 1e-6;     % F
v0 = 0;        % V
T  = 20e-3;    % s

%
vinfty = R*IS;
tau    = R*C;

%
t  = [0:T/1000:T];  % griglia temporale per la visualizzazione

% tensione ai capi del condensatore
v  = vinfty + (v0 - vinfty)*exp(-t/tau);

% corrente nel condensatore
iC = -C/tau*(v0 - vinfty)*exp(-t/tau);

% plot di v = v(t)
figure
plot(t, v, 'b-')
xlabel('t [s]')
ylabel('v(t) [V]');

% plot di i_C = i_C(t)
figure
plot(t, iC, 'r-')
xlabel('t [s]')
ylabel('i_C(t) [A]');
%
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Interpolazione
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% misure sperimentali
Ci = 1e-6*[0.1; 1.1; 0.5; 1.65; 2.01; 1.3; 0.9];
Qi = 1e-9*[3.8; 0.02; 0.2; 5.4; 3.2; 1.5; 2.8]; 
N  = length(Ci);

% riordinamento del vettore delle capacita' con il comando [X,I]=sort(Ci) 
[X,I] = sort(Ci);

% numero dei dati: N = n + 1, dove n e' il grado del polinomio interpolatore
N     = length(Ci);
n     = N - 1;

% riordinamento delle cariche
for k = 1:N
    Y(k,1) = Qi(I(k));
end

% griglia fine per la valutazione del polinomio interpolatore di grado n
xx = [X(1):(X(end)-X(1))/100:X(end)];

% interpolazione dei dati (carica misurata) con un polinomio di grado n
% (stabilizzazione con centering and scaling transformation (vedi help polyfit)
[pn, S, mu] = polyfit(X, Y, n);
pp          = polyval(pn, xx, [], mu);

% plot dei dati e del polinomio interpolante
figure
plot(X*1e6, Y*1e9, 'ro', xx*1e6, pp*1e9, 'b-')
xlabel('C [\mu F]')
ylabel('Q [nC]')

% calcolo della tensione a transitorio esaurito sulla base dei dati misurati (V = Q/C)
V = Y./X; 
 
% interpolazione dei dati (tensione determinata sulla base di misure sperimentali Q-C) con un polinomio di grado n
% (stabilizzazione con centering and scaling transformation (vedi help polyfit) 
[cn, S, mu] = polyfit(X, V, n);
pc          = polyval(cn, xx, [], mu);

% plot dei dati e del polinomio interpolante
figure
plot(X*1e6, V, 'ko', xx*1e6, pc, 'b-')
xlabel('C [\mu F]')
ylabel('V_\infty [Volt]')
%
return
