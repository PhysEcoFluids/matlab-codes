%% A script for the two island Lotka Volterra
close all
clear all
% Four cases are run at the same time with different forcing strengths
%% physical parameters
% define the time step to be a fraction of the forcing time
dt = 1e-4;
% model paramaters 1 difference in prey birth
% alpha1=0.75;beta1=4/3;mygam1=1;delta1=1;
% alpha2=2/3;beta2=4/3;mygam2=1.05;delta2=1;
% model paramaters 2 difference in hunting efficiency
alpha1 = 2/3; beta1 = 4/3; mygam1 = 1;delta1 = 1;
alpha2 = 2/3; beta2 = 4/3; mygam2 = 1;delta2 = 1/4;
% model paramaters 3 difference predator mortality
%  alpha1=2/3;beta1=4/3;mygam1=1;delta1=1;
%  alpha2=2/3;beta2=1/3;mygam2=0.1;delta2=1/4;

mydiff = [0 0.01 0.1 1]';
% population in thousands
%% paramaters for time stepping and analysis
numsteps = 200;
numouts = 1e4;

%% these are the functions for the DEs with prey moving
fprey1=@(t,prey1,pred1) alpha1*prey1-beta1*prey1.*pred1;
fpred1=@(t,prey1,pred1,pred2) -mygam1*pred1+delta1*prey1.*pred1-mydiff.*(pred1-pred2);
fprey2=@(t,prey2,pred2) alpha2*prey2-beta2*prey2.*pred2;
fpred2=@(t,prey2,pred2,pred1) -mygam2*pred2+delta2*prey2.*pred2-mydiff.*(pred2-pred1);

%% Initialization
preys1 = zeros(numouts+1,4);
preds1 = zeros(numouts+1,4);
preys2 = zeros(numouts+1,4);
preds2 = zeros(numouts+1,4);
ts = zeros(numouts+1,1);
prey1 = 20*ones(4,1); pred1 = 0.1*ones(4,1); t=0;
prey2 = prey1; pred2 = 0.1*ones(4,1); 
preys1(1,:) = prey1; preds1(1,:) = pred1; ts(1) = t;
preys2(1,:) = prey2; preds2(1,:) = pred2;
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
    preytilde1 = prey1+dt*fprey1(t,prey1,pred1);
    predtilde1 = pred1+dt*fpred1(t,prey1,pred1,pred2);
    preytilde2 = prey2+dt*fprey2(t,prey2,pred2);
    predtilde2 = pred2+dt*fpred2(t,prey2,pred2,pred1);
    ttilde = t+dt;
    prey1 = prey1+0.5*dt*(fprey1(t,prey1,pred1)+fprey1(ttilde,preytilde1,predtilde1));
    pred1 = pred1+0.5*dt*(fpred1(t,prey1,pred1,pred2)+fpred1(ttilde,preytilde1,predtilde1,predtilde2));
    prey2 = prey2+0.5*dt*(fprey2(t,prey2,pred2)+fprey2(ttilde,preytilde2,predtilde2));
    pred2 = pred2+0.5*dt*(fpred2(t,prey2,pred2,pred1)+fpred2(ttilde,preytilde2,predtilde2,predtilde1));
    t = ttilde;    
 end
 % store
 preys1(jj+1,:) = prey1; preds1(jj+1,:) = pred1; ts(jj+1) = t; 
 preys2(jj+1,:) = prey2; preds2(jj+1,:) = pred2;
end
%% Analysis and graphics
% define where to start the data so as to skip the transient
lhpt = floor(numouts/2);

% Figure 1 plots the phase portrait
figure(1)
clf
betterplots
for ii = 1:4
    subplot(2,2,ii)
    plot(preds1(lhpt:end,ii),preys1(lhpt:end,ii),'k-',preds2(lhpt:end,ii),preys2(lhpt:end,ii),'b-')
end
subplot(2,2,1)
ylabel('prey')
subplot(2,2,3)
ylabel('prey')
xlabel('pred')
subplot(2,2,4)
xlabel('pred')

figure(2)
clf
betterplots
% FFT based windowed spectra of prey
prey1now = preys1(lhpt:numouts,:);
prey2now = preys2(lhpt:numouts,:);
sz = size(prey1now);
dummy = floor(sz(1)/2);
mywin = ([1:dummy dummy:-1:0]')*[1 1 1 1]; % create 4 copies of a triangtle window
preys1f = fft(prey1now.*mywin,[],1);
preys2f = fft(prey2now.*mywin,[],1);
specs1 = abs(preys1f).^2;
specs2 = abs(preys2f).^2;
dom = 2*pi/(ts(end)-ts(lhpt))
% Set the number of frequencies to show and define those frequencies
numoms = 100;
oms = (0:numoms)*dom;
for ii = 1:4
    maxpsd(ii) = max(specs1(1:numoms+1,ii));
    subplot(2,2,ii)
    plot(oms,specs1(1:numoms+1,ii)/maxpsd(ii),'k-',oms,specs2(1:numoms+1,ii)/maxpsd(ii),'b-')
    title(['Max PSD = ' num2str(maxpsd(ii),4)])
end
subplot(2,2,1)
ylabel('scaled PSD')
subplot(2,2,3)
ylabel('scaled PSD')
xlabel('scaled frequency')
subplot(2,2,4)
xlabel('scaled frequency')

figure(3)
clf
betterplots
plot(ts(lhpt:end),preds1(lhpt:end,1)-preds2(lhpt:end,1),'r')
hold on
plot(ts(lhpt:end),preds1(lhpt:end,2)-preds2(lhpt:end,2),'g')
plot(ts(lhpt:end),preds1(lhpt:end,3)-preds2(lhpt:end,3),'m')
plot(ts(lhpt:end),preds1(lhpt:end,4)-preds2(lhpt:end,4),'c')
grid on
ylabel('pred difference')
xlabel('t')