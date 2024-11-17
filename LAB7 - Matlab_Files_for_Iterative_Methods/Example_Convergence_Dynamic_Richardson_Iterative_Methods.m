% This script studies the convergence of dynamic Richardson methods to solve a linear system

close all
clear all
clc

%
format
%

n     = 10;
B     = rand(10);

% gamma is the weight of the diagonal matrix added to A = B'*B
% select gamma = 0, gamma = 20, gamma = 40
gamma = 0;
A     = B'*B + gamma*eye(n);
% 

% check if A = A'
disp('check if A is symmetric')
A == A'
pause
clc

% check if A is spd
disp('check if A is positive definite')
eig(A) > 0
pause
clc

% solution of the system
xex   = ones(n,1);
b     = A*xex;

% Iteration parameters
toll        = 1e-12;
x0          = zeros(n,1);
nitmax      = 1000000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution of the linear system
% termination test: control of the increment (stop_test = 1)
%                   control of the residual  (stop_test = 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stop_test   = 1;

% choice of the preconditioner
P           = eye(n);

fprintf('   DYNAMIC RICHARDSON METHOD\n\n')

[xR, errR, kR] = richardson_dyn(A, b, x0, P, toll, nitmax, stop_test);

fprintf('# of iterations = %d\n', kR);
pause

format long e
[xex, xR]

pause
clc

format short e

fact = cond(inv(P)*A);

% verification of the reliability of the computed solution
if (stop_test == 1)
   true_error = norm(xex - xR);
   fprintf('termination test: control of the increment\n\n')
else
   true_error = norm(xex - xR)/norm(xex);
   fprintf('termination test: control of the residual\n\n')
end
%
fprintf('2-norm of the true absolute error      = %d\n', true_error)
est_error  = errR(end);
fprintf('2-norm of the estimated absolute error = %d\n', est_error)   
%
fprintf('error amplification factor             = %d\n', fact);
fprintf('2-norm of the estimated relative error * amplification factor = %d\n', fact*est_error)
%
return
