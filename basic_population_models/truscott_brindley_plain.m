%% A script for the Truscott and Brindley P and Z model
close all
clear all
% Four cases are run at the same time with different forcing strengths
%% physical parameters
% define the time step to be a fraction of the forcing time
dt = 1e-2;
% model paramaters
bigk = 108; r = 0.3; Rm = 0.7; alpha = 5.7; mu = 0.012; gamma = 0.05;
% population in thousands
P0s = [5; 10; 20; 50];
%% paramaters for time stepping and analysis
numsteps = 2e1;
numouts = 1000;

%% these are the functions for the DEs
fP=@(t,P,Z) r*P.*(1-P/bigk)-Rm*Z.*P.^2./(alpha^2+P.^2);
fZ=@(t,P,Z) gamma*Rm*Z.*P.^2./(alpha^2+P.^2)-mu*Z;

%% Initialization
Ps = zeros(numouts+1,4);
Zs = zeros(numouts+1,4);
ts = zeros(numouts+1,1);
P = P0s; Z = 5*ones(4,1); t = 0;
Ps(1,:) = P; Zs(1,:) = Z; ts(1) = t;
%% Loops
% Outer loop is over cycles to store at.
for jj = 1:numouts
% Inner loop is over individual time steps
 for ii = 1:numsteps;
  % Euler time stepping for comparison
%   P = P+dt*fP(t,P,Z);
%   Z = Z+dt*fZ(t,P,Z);
%   t = t+dt;
  % Heun time stepping as in notes
     Ptilde = P+dt*fP(t,P,Z);
     Ztilde = Z+dt*fZ(t,P,Z);
     ttilde = t+dt;
     P = P+0.5*dt*(fP(t,P,Z)+fP(ttilde,Ptilde,Ztilde));
     Z = Z+0.5*dt*(fZ(t,P,Z)+fZ(ttilde,Ptilde,Ztilde));
     t = ttilde;    
 end
 % store
 Ps(jj+1,:) = P; Zs(jj+1,:)=Z;ts(jj+1)=t; 
end
%% Analysis and graphics
% define where to start the data so as to skip the transient
lhpt = floor(numouts/2);
% Figure 1 plots the time series
figure(1)
clf
% This is my personal graphics improvement script
% write your own or comment out
betterplots
for ii = 1:4
    subplot(2,2,ii)
    plot(ts(1:end),Zs(1:end,ii),'b-',ts(1:end),Ps(1:end,ii),'k-')
    grid on
end
subplot(2,2,1)
ylabel('Z (blue) and P (black)')
subplot(2,2,3)
ylabel('Z (blue) and P (black)')
xlabel('t')
subplot(2,2,4)
xlabel('t')

% Figure 2 plots the phase portrait
figure(2)
clf
betterplots
for ii = 1:4
    subplot(2,2,ii)
    plot(Ps(1:end,ii),Zs(1:end,ii),'k-')
    grid on
end
subplot(2,2,1)
ylabel('P')
subplot(2,2,3)
ylabel('P')
xlabel('Z')
subplot(2,2,4)
xlabel('Z')


