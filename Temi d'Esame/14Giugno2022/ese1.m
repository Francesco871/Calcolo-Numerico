clc
close all
clear all

format long e;

%2)
O = zeros(3,3);
a = [29 , 2 , 1 ; 2 , 6 , 1 ; 1 , 1 , (1/5)];
A = [a , O , O ; O , a , O ; O , O , a];
n = 9;
xex = ones(n,1);
b = A*xex;

%1a. Metodo di Jacobi:
%check if there are zeros on the diagonal of A
diag(A); %NB ho definito io ha e so che non ci sono zeri => P invertibile
PJ    = diag(diag(A));
BJ    = eye(n) - inv(PJ)*A;

norm1J = norm(BJ,1);
norm2J = norm(BJ,2);
norm_infJ = norm(BJ,'inf');

%CONVERGENZA:
%check if A is diagonally dominant 
for i = 1:n
       s(i) = norm(A(i,:),1) - abs(A(i,i)); 
end
[diag(A), s']; % A NON è diagonalmente dominante per righe

% controllo il raggio spettrale:
rhoJ  = max(abs(eig(BJ))); % raggio > 1 => NON convergo

%1b. Metodo di Gauss-Seidel:
% check if A is symmetric
A == A';
%check if A is positive definite 
eig(A);
%Dalla verifica: A è sdp => GS converge

PGS    = tril(A); 
BGS    = eye(n) - inv(PGS)*A;

norm1GS = norm(BGS,1);
norm2GS = norm(BGS,2);
norm_infGS = norm(BGS,'inf');

rhoGS  = max(abs(eig(BGS))); % anche dalla verifica con rho convergo

%2. Uso richardson_stat
toll=1e-12;
nitmax=10000;
x0 = zeros(n,1);
stop_test=2;

[x,err,niter] = richardson_stat(A,b,PGS,1,toll,nitmax); %NB: QUESTA è LA VECCHIA FUNCTION
%[xGS,err,niter] = richardson_stat(A,b,x0,PGS,1,toll,nitmax,stop_test); NUOVA

est_rel_err = err(end);

%3. Errore relativo effettivo
est_rel_errC = (norm(xex-x,2))/(norm(xex,2));

K = cond(A); %fattore di condizionamento

[est_rel_err , est_rel_errC];

return