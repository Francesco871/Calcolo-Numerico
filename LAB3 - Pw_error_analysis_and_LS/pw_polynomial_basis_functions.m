%
close all
clear all
clc
%
a   = 0; 
b   = 1;
% number of elements of Tau_h
Mh  = 3;
% degree of the FE space
r   = 2;  
% uniform mesh size
H   = (b-a)/Mh; 
% vertices of Tau_h
xv  = [a:H:b];
% grid for the plot
Nsubdivisions = 20;
h   = H/Nsubdivisions;
x   = [a:h:b]';
Nx  = numel(x);
% loop over the elements of Tau_h
setfonts;
%
for i = 1:Mh
    % basis function set
    phi_i = zeros(Nx, r+1);
    % endpoints of K_i
    xi    = xv(i);
    xip1  = xv(i+1);
    % evaluation of the r+1 local basis function set
    phi_i_loc = eval_local_FE_basis_functions(xi, xip1, h, r);
    % plot of the basis function set associated with element K_i
    row_start = (i-1)*Nsubdivisions + 1;
    row_end   = row_start + Nsubdivisions;
    row_list  = [row_start:row_end];
    phi_i(row_list, :) = phi_i_loc;
    %
    plot(x, phi_i);
    xlabel('x')
    title('\phi_j|_{K_i}')
    pause
    close all
end
%
return
