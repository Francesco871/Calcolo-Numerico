% This script studies the convergence of stationary Richardson methods to solve a linear system

close all
clear all
clc

%
format
%

n     = 10;
B     = rand(n);

% gamma is the weight of the diagonal matrix added to A = B'*B
% select gamma = 0, gamma = 20, gamma = 40
gamma = 40;
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

% check if A is diagonally dominant by rows
disp('check if A is diagonally dominant by rows')
for i = 1:n
    difference(i) = abs(A(i,i)) - sum(abs(A(i, [1:i-1, i+1:n])));
end
%
difference, 
pause
clc

% solution of the system
xex   = ones(n,1);
b     = A*xex;

% Iteration parameters
alpha       = 1;
toll        = 1e-12;
x0          = zeros(n,1);
nitmax      = 100000;

%%%%%%%%%%%%%%%%
% Jacobi method
%%%%%%%%%%%%%%%%
PJ    = diag(diag(A));
BJ    = eye(n)-inv(PJ)*A;

fprintf('   JACOBI METHOD  \n\n')

% check if there is a norm of BJ such that ||B_J||<1
disp('check if there is a norm of B < 1')
%
[nBJ, index] = min([norm(BJ, 1), norm(BJ, 2), norm(BJ, 'inf')]);
if (nBJ < 1)
   %
   if (index == 1)
      disp('the 1-norm of B is the smallest')
   elseif (index == 2)
      disp('the 2-norm of B is the smallest')
   elseif (index == 3)
      disp('the inf-norm of B is the smallest')
   end
else
   disp('all the B-norms are >= 1')
end

pause
clc

% spectral radius of BJ
rhoJ = max(abs(eig(BJ)));
fprintf('rho(B_J)  = %d\n', rhoJ);
fprintf('||B_J||   = %d\n', nBJ);
pause
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution of the linear system
% termination test: control of the increment (stop_test = 1)
%                   control of the residual  (stop_test = 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stop_test      = 2;

% if nBJ < 1, then we can estimate the # of iterations needed to make the error <= toll
if ((stop_test == 1) & (nBJ < 1))
   %
   x1   = BJ*x0 + inv(PJ)*b;
   kmin = round(log(toll*(1-nBJ)/norm(x1-x0))/log(nBJ) - 1);
   fprintf('ESTIMATED # of iterations = %d\n', kmin);
   pause
end

if (rhoJ < 1)
   [xJ, errJ, kJ] = richardson_stat(A, b, x0, PJ, alpha, toll, nitmax, stop_test);

   fprintf('   ACTUAL # of iterations = %d\n', kJ);
   pause

   format long e
   [xex, xJ]

   pause
   clc

   format short e

   % verification of the reliability of the computed solution
   if (stop_test == 1)
      if (nBJ < 1)
         fact = nBJ/(1 - nBJ);
      else
         fact = rhoJ;
      end
   else
      fact = cond(inv(PJ)*A);
   end

   fprintf('error amplification factor = %d\n', fact);

   if (stop_test == 1)
      fprintf('termination test: control of the increment\n\n')
      true_error = norm(xex - xJ);
      fprintf('2-norm of the true absolute error      = %d\n', true_error)
      est_error  = errJ(end);
      fprintf('2-norm of the estimated absolute error = %d\n', est_error)   
      fprintf('2-norm of the estimated absolute error * amplification factor = %d\n', fact*est_error)
   else
      fprintf('termination test: control of the residual\n\n')
      true_error = norm(xex - xJ)/norm(xex);
      fprintf('2-norm of the true relative error      = %d\n', true_error)
      est_error  = errJ(end);
      fprintf('2-norm of the estimated relative error = %d\n', est_error)   
      fprintf('2-norm of the estimated relative error * amplification factor = %d\n', fact*est_error)
   end
   %
end
%
clc

%%%%%%%%%%%%%%%%%%%%%%%
% Gauss-Seidel method
%%%%%%%%%%%%%%%%%%%%%%%
PGS   = tril(A);
BGS   = eye(n)-inv(PGS)*A;

fprintf('   GAUSS-SEIDEL METHOD  \n\n')

% check if there is a norm of BGS such that ||B_GS||<1
disp('check if there is a norm of B < 1')
%
[nBGS, index] = min([norm(BGS, 1), norm(BGS, 2), norm(BGS, 'inf')]);
if (nBGS < 1)
   %
   if (index == 1)
      disp('the 1-norm of B is the smallest')
   elseif (index == 2)
      disp('the 2-norm of B is the smallest')
   elseif (index == 2)
      disp('the inf-norm of B is the smallest')
   end
else
   disp('all the B-norms are >= 1')
end

pause
clc

% spectral radius of BGS
rhoGS = max(abs(eig(BGS)));
fprintf('rho(B_GS)  = %d\n', rhoGS);
fprintf('||B_GS||   = %d\n', nBGS);
pause
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution of the linear system
% termination test: control of the increment (stop_test = 1)
%                   control of the residual  (stop_test = 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stop_test      = 2;

% if nBGS < 1, then we can estimate the # of iterations needed to make the error <= toll
if ((stop_test == 1) & (nBGS < 1))
   %
   x1   = BGS*x0 + inv(PGS)*b;
   kmin = round(log(toll*(1-nBGS)/norm(x1-x0))/log(nBGS) - 1);
   fprintf('ESTIMATED # of iterations = %d\n', kmin);
   pause
end

if (rhoGS < 1)
   [xGS, errGS, kGS] = richardson_stat(A, b, x0, PGS, alpha, toll, nitmax, stop_test);

   fprintf('   ACTUAL # of iterations = %d\n', kGS);
   pause

   format long e
   [xex, xGS]

   pause
   clc

   format short e

   % verification of the reliability of the computed solution
   if (stop_test == 1)
      if (nBGS < 1)
         fact = nBGS/(1 - nBGS);
      else
         fact = rhoGS;
      end
   else
      fact = cond(inv(PGS)*A);
   end

   fprintf('error amplification factor = %d\n', fact);

   if (stop_test == 1)
      fprintf('termination test: control of the increment\n\n')
      true_error = norm(xex - xGS);
      fprintf('2-norm of the true absolute error      = %d\n', true_error)
      est_error  = errGS(end);
      fprintf('2-norm of the estimated absolute error = %d\n', est_error)   
      fprintf('2-norm of the estimated absolute error * amplification factor = %d\n', fact*est_error)
   else
      fprintf('termination test: control of the residual\n\n')
      true_error = norm(xex - xGS)/norm(xex);
      fprintf('2-norm of the true relative error      = %d\n', true_error)
      est_error  = errGS(end);
      fprintf('2-norm of the estimated relative error = %d\n', est_error)   
      fprintf('2-norm of the estimated relative error * amplification factor = %d\n', fact*est_error)
   end
   %
end
%

return
