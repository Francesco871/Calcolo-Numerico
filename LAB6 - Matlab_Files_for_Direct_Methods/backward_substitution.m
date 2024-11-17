function x = backward_substitution(U, y)
% This function solves the upper triangular system U*x = y with the backward substitution method.
%
% x = backward_substitution(U, y)
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
   error('Matrix U is singular'); 
end
%
for j = n:-1:2
    y(j)     = y(j)/U(j,j); 
    y(1:j-1) = y(1:j-1)-y(j)*U(1:j-1,j);
end
%
y(1) = y(1)/U(1,1);
%
x    = y;
%
return
