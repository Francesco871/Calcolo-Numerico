% This script is an example of stability analysis in the numerical solution of linear algebraic systems
%
clear all
close all
clc
format 
%
n = 20;
B = rand(n);

% A is symmetric positive definite
A = B'*B;

% let us check it
A == A', 
pause
eig(A) > 0

% evaluate K_2(A)
K2 = cond(A)

% exact solution of the linear system 
x  = ones(n,1);
% right-hand side of the linear system
b  = A*x;

% flag to select the type of perturbation analysis
% flag = 1 -----> delta A = 0, delta b ~= 0
% flag = 2 -----> delta b = 0, delta A ~= 0
% flag = 3 -----> delta A ~= 0 & delta b ~= 0

flag = 1;

setfonts;

%
switch flag

case 1

     wA = 0;
     wb = 5;

case 2

     wA = 1;
     wb = 0;
     
case 3
 
     wA = 1e-4;
     wb = 5;
      
end

% define the perturbations

deltaA = wA*eye(n);

V      = [1:n]';
deltab = wb*ones(n,1).*sin(pi*V/n).*cos(pi*V/n);

plot(V, b, 'bo-', V, b+deltab, 'm*-') 
legend('b', 'b + \delta b')

% generate the perturbed input data

Ap = A + deltaA;
bp = b + deltab;

% solve the perturbed system

xp = Ap \ bp;

% compute the relative error

dx  = x - xp;
err = norm(dx)/norm(x);

% verify that deltaA satisfies the condition norm(inv(A))*norm(deltaA) < 1
condition = norm(inv(A))*norm(deltaA);

check     = condition < 1;

% compute the amplification factor
if (check > 0)
   
   amp_fact = K2 / (1 - K2*norm(deltaA)/norm(A))


   % compute the a priori estimate of the error

   err_est  = amp_fact*(norm(deltaA)/norm(A) + norm(deltab)/norm(b));

   % ratio between the estimated error and the actual error

   ratio    = err_est / err;

   err, err_est, ratio

else
   
   disp('Too large perturbation!!!!!') 
   
   condition
   
end
%
return
