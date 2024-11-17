% Authors: R.Sacco and G.Chiaravalli
% This script computes the derivative of the data on the displacement
% repeated over three periods obtained with the BCG measures and imported 
%in the form of .mat file
%The formula:
%
% y' = (y(t_0+h) - y(t_0))/h        (*)
%
% can be also generalized for each point with the command diff() that computes the
% difference in magnitude between each element of the input vector. 
% N.B.: before running the script use the graphical command set path
% to add to the Matlab path the folder ./matlab_library

clear all
close all 


%load the .mat file containing the bcg data 
load('D_BCG.mat')
m=numel(D_BCG);

ds_periodic=[D_BCG,D_BCG,D_BCG];
n=numel(ds_periodic);

t =linspace(0,2.4,n);
dt=t(2)-t(1);
s=m;
for i=1:m
    yprime(i)=(ds_periodic(i+1)-ds_periodic(i))/dt;
     s=s+1;
end

yprime=[yprime,yprime,yprime];
s=m;
for i=1:m
    ysecond(i)=(yprime(i+1)-yprime(i))/dt;
     s=s+1;
end
ysecond=[ysecond,ysecond,ysecond];

setfonts;
figure
plot(t,ds_periodic,'linewidth', 3)
xlabel('t (s)')
ylabel('displacement (g cm)')

figure
plot(t,yprime,'linewidth', 3)
xlabel('t (s)')
ylabel('velocity (g cm s^-1)')

figure
plot(t,ysecond,'linewidth', 3)
xlabel('t (s)')
ylabel('acceleration (dyne)')






