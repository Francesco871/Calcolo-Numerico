function Q = quadrature(xb, H, f, wn, xn)
% This function computes the numerical approximation of the integral of f=f(x) over the interval [a,b]
%
% Q = quadrature(xb, H, f, wn, xn)
%
M = numel(H);
Q = 0;
%
for i = 1:M
    h  = H(i);
    fn = f(xb(i) + h/2*xn);
    Q  = Q + h/2*sum(wn.*fn);
end
%
return
