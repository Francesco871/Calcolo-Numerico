function [a, B] = least_squares(x, y, n)
% This function computes the coefficients of the polynomial p_n^* that approximates the data pairs (x_i, y_i),
% i=0, 1, ..., m, in the least squares sense
%
% [a, B] = least_squares(x, y, n)
%
B = zeros(n+1);
b = zeros(n+1,1);
%
for i = 0 : n,
    for j = 0 : n
        B(i+1, j+1) = (x.^i)'*x.^j;
    end
    b(i+1) = y'*x.^i;
end
%
a = B \ b;
%
return
