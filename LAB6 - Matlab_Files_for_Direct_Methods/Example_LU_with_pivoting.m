% This script: 
% 1. computes the matrices P, L and U that verify P*A = L*U using user-coded functions and the L, U, P Matlab command. 
% 2. computes the residual matrix R = P*A - L*U that ideally should be equal to the zero matrix.
% 3. solves the linear system A*x = b using the user-coded fwd and bckwd substitution algorithms and the \ Matlab command
%
clear all
close all
clc
%
format short 

% 
A = [1 2 3; 2 4 5; 7 8 9];
n = max(size(A));

% Construct the rhs corresponding to the solution x = [1, 1, 1]^T
xex = ones(n,1);
b   = A*xex;

% Check that A is not singular
d = det(A);
fprintf('det(A) = %d\n', d) 
pause

% Check the necessary and sufficient condition for A to admit a unique LU factorization with L_{ii} = 1, i=1,...,n
for i = 1:(n-1)
    det(A(1:i, 1:i)), pause
end

% By using the switch ifact we can decide whether to run the case 1 (user coded) or 2 (matlab function)
ifact = 1;

switch ifact
    case 1
%
[L1, U1, P1] = lu_factorization_with_pivoting(A)
R1           = P1*A - L1*U1,
pause

    case 2
% Factorization with Matlab built-in function
[L1, U1, P1] = lu(A)
R2           = P1*A - L1*U1,
pause
end
 

AP = P1*A;
for i = 1:(n-1)
    AP(1:i, 1:i), det(AP(1:i, 1:i)), pause
end

% Solve the two triangular linear systems using the user-coded substitution algorithms

% Solve the linear system to compute x
% y = L\b;
% x = U\yx = backward_substitution(U, y)

% forward substitution
y           = forward_substitution(L1, P1*b);
% backward substitution
x           = backward_substitution(U1, y)
% Solve the linear system with the \ command = LU factorization with pivoting + fwd + bckwd substitutiom
xMatlab     = A \ b
%
return
