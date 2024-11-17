clc
close all
clear all

n = 5;
A = pascal(n);

xex = ones(n,1);
b = A*xex;

%Gauss-Seidel
% check if A is symmetric
A == A';
%check if A is positive definite 
eig(A);
%Dalla verifica: A è sdp => GS converge

PGS    = tril(A); 
BGS    = eye(n) - inv(PGS)*A;

toll=1e-12;
nitmax=100000;
x0 = zeros(n,1);
stop_test=2;

rhoGS     = max(abs(eig(BGS)))

kiter     = log(toll)/log(rhoGS);

[xGS,errGS,niterGS] = richardson_stat(A,b,PGS,1,toll,nitmax); %NB: QUESTA è LA VECCHIA FUNCTION
%[xGS,err,niter] = richardson_stat(A,b,x0,PGS,1,toll,nitmax,stop_test); NUOVA

true_rel_Err = (norm(xex-xGS,2))/(norm(xex,2));

k=cond(A);

return