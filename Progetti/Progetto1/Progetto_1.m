clc
close all
clear all

setfonts;
format short e;

%2)IMPLEMENTAZIONE E VISUALIZZAZIONE

%1.Assegno i valori
v0 = 0;
Is = 1e-06;
R = 2e03;
C = 1e-06;
T = 20e-03;

%2.Definisco colonna t
t = [0:(T/1000):T]';

%3.Calcolo v(t) in ogni istante t
v_inf = R*Is;
tau = R*C;
v = v_inf+(v0-v_inf)*exp(-t/tau);

%4.Plotto v
figure
plot(t, v, 'b-')
xlabel('t[s]')
ylabel('v(t)[V]')

%5.Calcolo i_c in ogni istante t
i_c = -(C/tau)*(v0-v_inf)*exp(-t/tau);

%6.Plotto i_c
figure
plot(t, i_c, 'r-')
xlabel('t[s]')
ylabel('i_c(t)[A]')

%3)TRATTAMENTO DEI DATI

%1.Costruisco vettori C_i e Q_i
C_i = [0.1; 1.1; 0.5; 1.65; 2.01; 1.3; 0.9]*1e-06;
Q_i = [3.8; 0.02; 0.2; 5.4; 3.2; 1.5; 2.8;]*1e-09;

%2.Utilizzo "sort" per ordinare C_i in modo crescente in X e salvo indici
%in I
[X,I] = sort(C_i);

%3.Costruisco Y con "for" e I
N = numel(C_i);
for i=1:N
    Y(i,1)=Q_i(I(i));
end

%4.Definisco xx tra X(1) e X(end) con spaziatura uniforme
xx = [X(1):(X(end)-X(1))/100:X(N)];

%5.Calcolo i coefficienti del polinomio che interpola Y in X con grado n=N-1
n = N-1;
cn = polyfit(X,Y,n);
c0 = fliplr(cn);

%6.Valuto il polinomio pp in ogni punto di xx e poi lo plotto
pp = polyval(cn, xx);

figure
plot(xx*1e06, pp*1e09, 'b-', X*1e06, Y*1e09, 'ok')
xlabel('C[\muF]')
ylabel('Q[nC]')

%7.Calcolo v_inf per ogni coppia X,Y e lo salvo in V
V = Y./X;

%8.Calcolo i coefficienti del polinomio che interpola V in X
cn_tilde = polyfit(X,V,n);
c0_tilde = fliplr(cn_tilde);

%9.Valuto il polinomio pc in ogni punto di xx e poi lo plotto
pc = polyval(cn_tilde, xx);

figure
plot(xx*1e06, pc, 'b-', X*1e06, V, 'ok')
xlabel('C[\muF]')
ylabel('v_\infty[V]')

return