% This script is an example of the use of LU factorization to solve a system where
% A is singular but such that it admits a unique LU factorization
% In this case, A is a singular matrix but verifies the existence and uniqueness theorem for the LU factorization
clear all
close all 
clc
%
format
%
n        = 10;
A        = rand(n);
A        = A'*A;
A(end,:) = 0;
%
d        = det(A);
fprintf('det(A) = %d\n', d) 
pause

% check the necessary and sufficicnt condition for A to admit a unique LU factorization with L_{ii}=1,...,n
for k = 1:(n-1)
    d = det(A(1:k, 1:k))
    pause
end

% characterize Im(A)
w = rand(n,1);
%
A*w

% Factorization with a user-coded function
[L, U] = lu_factorization(A)

% Check that A == LU 
R = A - L*U 

% Construct the right-hand side b in such a way that the exact solution is x = [1; 1; 1;, .....; 1 ; alpha]
xex    = ones(n,1);
alpha  = 10;
xex(n) = alpha;
b      = A*xex;

% Check that b \in Im(A)
b(n) == 0

% Solve the linear system to compute x
% L*y = b;
% U*x = y

% Forward substitution
y   = forward_substitution(L, b);

% Backward substitution
x   = backward_substitution_singular(U, y, alpha);

% Show the solution
format long e
[xex, x]

return
