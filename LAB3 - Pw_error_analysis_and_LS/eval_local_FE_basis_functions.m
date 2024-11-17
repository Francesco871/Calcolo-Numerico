function phi_i_loc = eval_local_FE_basis_functions(xn, x, r)
% This function computes the local basis functions of degree r, r >=1,
% over the element K_i = [x_i, x_{i+1}] belonging to the triangulation Tau_h
%
% phi_i_loc = eval_local_FE_basis_functions(xn, x, r)
%
% construction of the basis functions of degree r>=1 over K_i
for j = 1:r+1
    % construction of the local basis function phi_j(x). j = 1, ..., r+1
    pol = 1;
    xj  = xn(j);
    %
    for k = 1:r+1
        xk = xn(k);
        %
	if (k ~= j)
            pol  = pol.*(x - xk)/(xj - xk);
	end
    end
    phi_i_loc(:,j) = pol;
end
% 
return
