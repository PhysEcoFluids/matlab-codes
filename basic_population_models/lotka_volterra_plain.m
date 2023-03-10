%% This code solves the standard Lotka Volterra model
close all
clear all
% Four cases are run at the same time with different forcing strengths
%% physical parameters
% define the time step to be a fraction of the forcing time
dt = 1e-3;
% model paramaters
alpha = 2/3;beta = 4/3;mygam = 1;delta = 1;
% population in thousands
prey0s = [0.1; 0.5; 5; 10];
%% paramaters for time stepping and analysis
numsteps = 1e1;
numouts = 10000;

%% these are the functions for the DEs
fprey=@(t,prey,pred) alpha*prey-beta*prey.*pred;
fpred=@(t,prey,pred) -mygam*pred+delta*prey.*pred;

%% Initialization
preys = zeros(numouts+1,4);
preds = zeros(numouts+1,4);
ts = zeros(numouts+1,1);
prey = prey0s; pred = 0.1*ones(4,1); t = 0;
preys(1,:) = prey; predss(1,:) = pred; ts(1) = t;
%% Loops
% Outer loop is over cycles to store at.
for jj=1:numouts
% Inner loop is over individual time steps
 for ii=1:numsteps;
  % Euler time stepping for comparison uncomment the part below
%   x = x+dt*fx(t,x,v);
%   v = v+dt*fv(t,x,v);
%   t = t+dt;
  % Heun time stepping
    preytilde = prey+dt*fprey(t,prey,pred);
    predtilde = pred+dt*fpred(t,prey,pred);
    ttilde = t+dt;
    prey = prey+0.5*dt*(fprey(t,prey,pred)+fprey(ttilde,preytilde,predtilde));
    pred = pred+0.5*dt*(fpred(t,prey,pred)+fpred(ttilde,preytilde,predtilde));
    t = ttilde;    
 end
 % store
 preys(jj+1,:) = prey; preds(jj+1,:) = pred;ts(jj+1) = t; 
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
    plot(ts(lhpt:end),preds(lhpt:end,ii),'b-',ts(lhpt:end),preys(lhpt:end,ii),'k-')
    grid on
end
subplot(2,2,1)
ylabel('pred and prey (black)')
subplot(2,2,3)
ylabel('pred and prey (black)')
xlabel('t')
subplot(2,2,4)
xlabel('t')

% Figure 2 plots the phase portrait
figure(2)
clf
betterplots
for ii = 1:4
    subplot(2,2,ii)
    plot(preds(lhpt:end,ii),preys(lhpt:end,ii),'k-')
    grid on
end
subplot(2,2,1)
ylabel('prey')
subplot(2,2,3)
ylabel('prey')
xlabel('pred')
subplot(2,2,4)
xlabel('pred')


