% Questo script contiene la soluzione del progetto #2
clear all
close all
format short e
%
setfonts;
% valori dei parametri circuitali
V_in = 1;
R_in = 100;
R_1  = R_in;
R_2  = 2*R_in;
R_3  = 3*R_in;
R_4  = 4*R_in;
%
G_in = 1/R_in;
G_1  = 1/R_1;
G_2  = 1/R_2;
G_3  = 1/R_3;
G_4  = 1/R_4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Analisi circuitale con metodi diretti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y    = zeros(6);
b    = zeros(6,1);

% matrice delle ammettenze
Y(1,2) = -G_1;
Y(1,4) = -G_2;
Y(1,1) = -sum(Y(1,:)) + G_in;
%
Y(2,1) = -G_1;
Y(2,3) = -G_1;
Y(2,5) = -G_2;
Y(2,2) = -sum(Y(2,:));
%
Y(3,2) = -G_1;
Y(3,6) = -G_2;
Y(3,3) = -sum(Y(3,:));
%
Y(4,1) = -G_2;
Y(4,5) = -G_3;
Y(4,4) = -sum(Y(4,:)) + G_4;
%
Y(5,2) = -G_2;
Y(5,4) = -G_3;
Y(5,6) = -G_3;
Y(5,5) = -sum(Y(5,:)) + G_4;
%
Y(6,3) = -G_2;
Y(6,5) = -G_3;
Y(6,6) = -sum(Y(6,:)) + G_4; 

% verifica esistenza e unicita' della fattorizzazione LU senza pivoting
for k = 1:5
    det(Y(1:k, 1:k))
end

% fattorizzazione LU senza pivoting
[L, U] = lu_factorization(Y)
pause

% risoluzione sistema lineare

% termine noto
b(1)   = V_in*G_in;

% risoluzione sistema triangolare inferiore
y      = forward_substitution(L, b);
% risoluzione sistema triangolare superiore
xc     = backward_substitution(U, y);
% confronto con il comando \
xm     = Y \ b;
%
format long e
[xm, xc]

% valutazione differenza delle due soluzioni
format short e
norm(xm - xc, 'inf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Analisi circuitale con metodi iterativi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verifica Y simmetrica e definita positiva
Y == Y'
pause
eig(Y)
pause
%

% risoluzione sistema lineare con ilo metodo di Gauss-Seidel
x0        = zeros(6,1);
alpha     = 1;
P         = tril(Y);
toll      = 1e-12;
nitmax    = 100000;
stop_test = 2; 

[xGS, errGS, kGS] = richardson_stat(Y, b, x0, P, alpha, toll, nitmax, stop_test);

kGS, 

format long e
[xm, xGS]

format short e

% errore relativo effettivamente commesso
true_err_rel = norm(xm - xGS)/norm(xm)

% errore relativo stimato
est_err_rel  = errGS(end)

% fattore di condizionamento
kappa        = cond(inv(P)*Y)

return
