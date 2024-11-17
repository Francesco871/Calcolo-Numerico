function [l_i, x] = eval_global_l_i(a, b, n, xn)
% [l_i, x] = eval_l_i(a, b, n, xn)
h   = (b-a)/100;
x   = [a:h:b]';
N   = numel(x);
l_i = zeros(N, n+1);
%
for i = 1:n+1
    pol = 1;
    xi  = xn(i);
    for j = 1:n+1
        xj = xn(j);
        %
        if (j~=i)
           pol = pol.*(x - xj)/(xi - xj);
        end
    end
    l_i(:,i) = pol;
end
%
return
