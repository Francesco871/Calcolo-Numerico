clc
close all
clear all

format short e;

a = 0;
b = 2*pi;
f = @(x) x.*exp(-0.1.*x.^2);
Iex = 5*(1-exp(-(2*pi^2)/5));
Mh = [1, 2, 4, 8, 16, 32, 64, 128, 256, 528];
n = numel(Mh);
for i = 1:n
    h  = (b-a)/Mh(i);
    xv = [a:h:b];
    xb = (xv(2:end) + xv(1:end-1))/2;
    H  = diff(xv);

%1. Errore
%Punto Medio
    wnMP = 2;
    xnMP = 0;
    QMP(i)  = quadrature(xb, H, f, wnMP, xnMP);
    errMP(i) = abs(Iex - QMP(i));

%Trapezio
    wnTR = [1;  1];
    xnTR = [-1; +1];
    QTR(i)  = quadrature(xb, H, f, wnTR, xnTR);
    errTR(i) = abs(Iex - QTR(i));

%Cavalieri-Simpson
    wnCS = 2*[1/6; 4/6; 1/6];
    xnCS = [-1; 0; +1];
    QCS(i)  = quadrature(xb, H, f, wnCS, xnCS);
    errCS(i) = abs(Iex - QCS(i));

%Gauss-Legendre
    wnGL = [1; 1];
    xnGL = [-1/sqrt(3); +1/sqrt(3)];
    QGL(i)  = quadrature(xb, H, f, wnGL, xnGL);
    errGL(i) = abs(Iex - QGL(i));

end

%2. Ordine di Convergenza
    h_p = 1./Mh;
    %p = (log(err(end-1))-log(err(end)))./(log(h_p(end-1))-log(h_p(end)));
    pMP = (log(errMP(end-1))-log(errMP(end)))./(log(h_p(end-1))-log(h_p(end)));
    pTR = (log(errTR(end-1))-log(errTR(end)))./(log(h_p(end-1))-log(h_p(end)));
    pCS = (log(errCS(end-1))-log(errCS(end)))./(log(h_p(end-1))-log(h_p(end)));
    pGL = (log(errGL(end-1))-log(errGL(end)))./(log(h_p(end-1))-log(h_p(end)));

return