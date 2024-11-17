function y = forward_substitution(L, b)
% This function solves the lower triangular system L*y = b with the forward substitution method.
%
% y = forward_substitution(L, b)
%
[n,m] = size(L);
%
if (n ~= m)
   clc
   error('Only square systems'); 
end
%
if (min(abs(diag(L))) == 0)
   clc
   error('Matrix L is singular'); 
end
%
for j = 1:n-1
    b(j)     = b(j)/L(j,j); 
    b(j+1:n) = b(j+1:n)-b(j)*L(j+1:n,j);
end
%
b(n) = b(n)/L(n,n);
%
y    = b;
%
return
