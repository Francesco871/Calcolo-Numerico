% This script computes the matrices L and U that verify A = L*U using a user-coded function
% and solves the linear system using forward and backward substitution with user-coded functions
%
clear all
close all
clc
% this command resets the output format on the screen
format 
% 
flag = 4;

switch flag

case 1
% In this case, A verifies the existence and uniqueness theorem for the LU factorization
A = [2 1 3; 5 4 6; 7 8 0];

case 2
% In this case, A does not verify the existence and uniqueness theorem for the LU factorization 
% although A is nonsingular (so that the linear system is uniquely solvable)
A = [1 2 3; 2 4 5; 7 8 9];

case 3
% In this case, A is a symmetric positive definite (sdp) matrix
A = rand(6);
A = A'*A;
eig(A), pause

case 4
% In this case, A is a singular matrix but verifies the existence and uniqueness theorem for the LU factorization
A        = rand(6);
A        = A'*A;
A(end,:) = 0;
A(:,end) = 0;
%
d        = det(A);
fprintf('det(A) = %d\n', d) 
pause

end
%
n = max(size(A));

% Check the necessary and sufficient condition for A to admit a unique LU factorization with L_{ii} = 1, i=1,...,n
for i = 1:(n-1)
    det(A(1:i, 1:i)), pause
end

% Factorization with a user-coded function
[L, U] = lu_factorization(A)

% Check that A == LU 
R = A - L*U 

% Construct the right-hand side b in such a way that the exact solution is x = [1; 1; 1]
xex = ones(n,1);
b   = A*xex;

% Solve the linear system to compute x
% L*y = b;
% U*x = y

% Check whether A is nonsingular
det(A), pause

% Forward substitution
y   = forward_substitution(L, b);

% Backward substitution
x   = backward_substitution(U, y);

% Show the solution
format long e
x
%
return
