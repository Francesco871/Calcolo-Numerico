clc
close all
clear all

setfonts;

%2)ANALISI CIRCUITALE CON METODI DIRETTI
Vin = 1; %[V]
Rin = 100; %[ohm]

R1 = 100;
R2 = 200;
R3 = 300;
R4 = 400;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;

%1.Introduco e visualizzo la matrice delle ammettenze Y e il termine noto b
Y = [(G1+G2+1/Rin) , (-G1) , 0 , (-G2), 0 , 0 ;
    (-G1) , (2*G1+G2) , (-G1) , 0 , (-G2) , 0 ;
    0 , (-G1) , (G1+G2) , 0 , 0 , (-G2) ;
    (-G2) , 0 , 0 , (G2+G3+G4) , (-G3) , 0 ;
    0 , (-G2) , 0 , (-G3) , (G2+2*G3+G4) , (-G3);
    0 , 0 , (-G2) , 0 , (-G3) , (G2+G3+G4)];

b = [(Vin/Rin) ; 0 ; 0 ; 0 ; 0 ; 0];

%2.Verifico l'esistenza ed unicità della fattorizzazione LU di Y
%NBverificata se e solo se il det dei minori principali di Y è diverso da 0
format short e;

n = max(size(Y));
for i = 1:(n-1)
    det_minors(i) = det(Y(1:i, 1:i));
end

%3.Calcolo la fattorizzazione LU utilizzando la function lu_factorization
[L, U] = lu_factorization(Y);

%4.Si consideri la soluzione del problema Yx=b
format long e;

%4a.Utilizzo la fattorizzazione LU di Y per risolvere il sistema con le funzioni 
%   forward_substitution e backward_substitution, memorizzo in xc la soluzione 
%   del sistema triangolare superiore
yc = forward_substitution(L, b);
xc = backward_substitution(U, yc);

%4b.Risolvo il sistema utilizzando il comando \ di Matlab, memorizzo in xm il risultato
ym = L\b;
xm = U\ym;

%5.Calcolo la norma infinito della differenza tra xm e xc 
format short e;

diff = norm(xm-xc, 'inf');

%3)ANALISI CIRCUITALE CON METODI ITERATIVI

%1.Verifico che la matrice Y è simmetrica e definita positiva
Y == Y';
eig(Y);

%2.Utilizzo la function richardson_stat per risolvere il sistema con il metodo di
%  Gauss-Seidel e memorizzo la soluzione calcolata nel vettore xGS
format long e;
x0 = zeros(6,1);
toll = 1e-12;
nitmax = 10000;
stop_test = 2;

PGS = tril(Y);
[xGS,err,k] = richardson_stat(Y,b,x0,PGS,1,toll,nitmax,stop_test);

%3.
format short e;

%3a.Calcolo l'errore relativo effettivamente commesso
true_rel_err = (norm(xm-xGS,2))/(norm(xm,2));

%3b.Considero l'errore relativo stimato est_rel_err dato dall'ultima componente del vettore errore 
%   restituito in uscita dalla function richardson_stat
est_rel_err = err(k);

%4.Calcolo il fattore di condizionamento dell’errore relativo
fact = cond(inv(PGS)*Y);

return