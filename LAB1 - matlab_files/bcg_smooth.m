% Authors: R.Sacco and G.Chiaravalli
% This script computes the derivative of the data on the displacement
% obtained with the BCG measures and imported in the form of .mat file
%The formula:
%
% y' = (y(t_0+h) - y(t_0))/h        (*)
%
% can be also generalized for each point with the command diff() that computes the
% difference in magnitude between each element of the input vector. 
% The command smooth is also applied to remove the discontinuities in the curves 
% N.B.: before running the script use the graphical command set path
% to add to the Matlab path the folder ./matlab_library


clear all
close all 

%load the .mat file containing the bcg data 
load('D_BCG.mat')
n=numel(D_BCG);
%the data on the displacement are referred to a time interval from 0s - 0.8s
t       =linspace(0,0.8,n);
t_der   =linspace(0,0.8,n-1);

%VELOCITY-first derivative -----------------
for i=1:n-1
    yprime(i)=(D_BCG(i+1)-D_BCG(i))/(t(i+1)-t(i));
end
%interpolation on t nodes
yprime=interp1(t_der,yprime,t);
% Smoothing for velocity
V_BCG_smooth = smooth(yprime,100,'lowess')';

%ACCELERATION-second derivative -----------------
for i=1:n-1
    ysecond(i)=(yprime(i+1)-yprime(i))/(t(i+1)-t(i));
end
%interpolation on t nodes
ysecond=interp1(t_der,ysecond,t);
% Smoothing for acceleration
A_BCG_smooth = smooth(ysecond,100,'lowess')';

setfonts;
%plot pf displacement
figure
plot(t,D_BCG,'linewidth',3)
xlabel('t(s)')
ylabel('displacement(g cm)')

%plot pf velocity
figure
plot(t,yprime,'linewidth',3), hold on
plot(t,V_BCG_smooth,'r--','linewidth',3)
xlabel('t(s)')
ylabel('velocity (g cm s^-1)')

%plot of acceleration
figure
plot(t,ysecond,'linewidth',3), hold on
plot(t,A_BCG_smooth,'r--','linewidth',3)
xlabel('t (s)')
ylabel('acceleration (dyne)')


%plot of acceleration
figure
plot(t,ysecond,'linewidth',3), hold on
plot(t,A_BCG_smooth,'r--','linewidth',3)
xlabel('t(s)')
ylabel('acceleration (dyne)')
axis([0 0.8 -70000 70000])

