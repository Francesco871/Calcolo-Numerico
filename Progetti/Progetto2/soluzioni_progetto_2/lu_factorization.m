function [L, U] = lu_factorization(A)
%
% This function computes the LU factorization of a given matrix A.
% This process consists of determining two triangular matrices L and U such that:
%
% 1. L*U = A
% 2. L_{ii} = 1, i=1, ..., n
%
% [L, U] = lu_factorization(A)
%
[n,m] = size(A);
%
if (n~= m) 
   error('Only square systems'); 
end
%
for k = 1:n-1
    if (A(k,k) == 0)
       clc
       fprintf('\n')
       fprintf('Factorization step k = %d\n',k)
       fprintf('Pivot                = %d\n',A(k,k))
       error('Null pivot element'); 
    end
    %
    A(k+1:n,k) = A(k+1:n,k)/A(k,k);
    %
    for j = k+1:n
        i      = [k+1:n]; 
        A(i,j) = A(i,j)-A(i,k)*A(k,j);
    end
end
%
L = tril(A,-1) + eye(n);
U = triu(A);
%
return
