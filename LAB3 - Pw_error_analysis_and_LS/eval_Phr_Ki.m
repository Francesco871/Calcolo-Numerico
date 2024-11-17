function [xx, Phr_Ki, xn, fn] = eval_Phr_Ki(xi, xip1, r, fun)
% This function computes the interpolation polynomial of degree r, with r >=0, 
% of the restriction of a given function f=f(x) over an element 
% K_i = [xi, xip1] belonging to the triangulation Tau_h, i = 1, ..., M_h 
h  = xip1 - xi;
dh = h;
if (r==0)
% the FE space P0(K_i)
   xn = (xi+xip1)/2;
   fn = fun(xn);
   % grid for interpolation and plot
   xx = [xi:dh/20:xip1];  
   Phr_Ki = fn*ones(size(xx));	 
elseif (r > 0)
% the FE space Pr(K_i), r>=1
   dh = h/r;
   % finite element nodes on K_i
   xn = [xi:dh:xip1]; 
   % grid for interpolation and plot
   xx = [xi:dh/20:xip1];  
   % nodal values of f(x) over K_i
   fn = fun(xn);           
   % construction of the interpolating polynomial of degree r>=1 over K_i
   Phr_Ki = 0;
   for j = 1:r+1
   % construction of the local basis function phi_j(x). j = 1, ..., r+1
       phi_j = 1;
       xj    = xn(j);
       %
       for k = 1:r+1
           xk = xn(k);
           %
	   if (k ~= j)
               phi_j = phi_j.*(xx - xk)/(xj - xk);
	   end
       end
       % interpolation
       Phr_Ki = Phr_Ki + fn(j)*phi_j; 
   end
end
% 
return
