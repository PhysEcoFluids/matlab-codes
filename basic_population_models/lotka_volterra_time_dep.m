%% A script for the Lotka Volterra with a time varying parameter
close all
clear all
% Four cases are run at the same time with different forcing strengths
%% physical parameters
dt = 1e-3;
% model paramaters
alpha = 2/3; beta = 4/3; mygam = 1; delta = 1;
% time dependence parameters
varper = 20;
varamp = 0.4;
%varamp = 0.5;

% population in thousands
prey0s = [0.1; 0.5; 5; 10];
%% paramaters for time stepping and analysis
numsteps = 2e1;
numouts = 80000;

%% these are the functions for the DEs
fprey=@(t,prey,pred) alpha*(1+varamp*sin(2*pi*t/varper))*prey-beta*prey.*pred;
fpred=@(t,prey,pred) -mygam*pred+delta*prey.*pred;

%% Initialization
preys = zeros(numouts+1,4);
preds = zeros(numouts+1,4);
ts = zeros(numouts+1,1);
prey = prey0s; pred = 0.1*ones(4,1); t = 0;
preys(1,:) = prey; preds(1,:) = pred; ts(1) = t;
%% Loops
% Outer loop is over cycles to store at.
for jj = 1:numouts
% Inner loop is over individual time steps
 for ii = 1:numsteps;
  % Euler time stepping for comparison
%   x = x+dt*fx(t,x,v);
%   v = v+dt*fv(t,x,v);
%   t = t+dt;
  % Heun time stepping as in notes
    preytilde = prey+dt*fprey(t,prey,pred);
    predtilde = pred+dt*fpred(t,prey,pred);
    ttilde = t+dt;
    prey = prey+0.5*dt*(fprey(t,prey,pred)+fprey(ttilde,preytilde,predtilde));
    pred = pred+0.5*dt*(fpred(t,prey,pred)+fpred(ttilde,preytilde,predtilde));
    t = ttilde;    
 end
 % store
 preys(jj+1,:) = prey; preds(jj+1,:) = pred; ts(jj+1) = t; 
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
end
subplot(2,2,1)
ylabel('pred and prey')
subplot(2,2,3)
ylabel('pred and prey')
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
end
subplot(2,2,1)
ylabel('prey')
subplot(2,2,3)
ylabel('prey')
xlabel('pred')
subplot(2,2,4)
xlabel('pred')

figure(4)
clf
betterplots
% FFT based windowed spectra of prey
preynow = preys(lhpt:numouts,:);
sz = size(preynow);
dummy = floor(sz(1)/2);
mywin = ([1:dummy dummy:-1:0]')*[1 1 1 1]; % create 4 copies of a triangular window
preysf = fft(preynow.*mywin,[],1);
specs1 = abs(preysf).^2;
dom = 2*pi/(ts(end)-ts(lhpt))
% Set the number of frequencies to show and define those frequencies
numoms = 100;
oms = (0:numoms)*dom;
for ii = 1:4
    maxpsd(ii) = max(specs1(1:numoms+1,ii));
    subplot(2,2,ii)
    plot(oms*varper/(2*pi),specs1(1:numoms+1,ii)/maxpsd(1),'k-')
    %title(['Max PSD = ' num2str(maxpsd(ii),4)])
    axis([0 2.5 0 0.5])
    grid on
end
subplot(2,2,1)
ylabel('scaled PSD')
subplot(2,2,3)
ylabel('scaled PSD')
xlabel('scaled frequency')
subplot(2,2,4)
xlabel('scaled frequency')
