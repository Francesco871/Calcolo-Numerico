function x = backward_substitution_singular(U, y, alpha)
% This function solves the upper triangular system U*x = y with the backward substitution method
% in the case where U(n, n) = 0
%
% x = backward_substitution_singular(U, y, alpha)
%
[m, n] = size(U);
%
if (n ~= m)
   clc
   error('Only square systems'); 
end
%
if (min(abs(diag(U))) == 0)
   clc
   disp('Matrix U is singular'); 
   pause
end
%
for j = n:-1:2
    if (U(j,j) == 0)
       y(j) = alpha;
    else
       y(j) = y(j)/U(j,j); 
    end
    y(1:j-1) = y(1:j-1)-y(j)*U(1:j-1,j);
end
%
y(1) = y(1)/U(1,1);
%
x    = y;
%
return
