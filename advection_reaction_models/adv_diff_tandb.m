%% This script solves the Turcott and Brindley model with advection and diffusion
clear all
close all
%% Definitions and parameters
% model paramaters
bigk = 108; r = 0.3; Rm = 0.7; alpha = 5.7; mu = 0.012; gamma = 0.05;
%define a spatial grid
xmin = -1e1;
xmax = 1e1;
N = 2^9;
x = linspace(xmin,xmax,N+1); x=x(1:end-1);
dx = x(2)-x(1);
% switches for plots and saving
LINE_PLOT = 1;
SAVE_FIELDS = 1;
%make initial condition
P = 15+20*(sech(((x+5)/1).^2)+sech(((x-5)/1).^2).*cos(2*pi*x/0.4));
%P = 5+30*(sech(((x+5)/1).^2)+sech(((x-5)/1).^2).*cos(2*pi*x/0.4));
Z = 5+2*rand(size(P));
P0 = P;
Z0 = Z;
%% these are the functions for the DEs
fP=@(t,P,Z) (r*P.*(1-P/bigk)-Rm*Z.*P.^2./(alpha^2+P.^2)).*(P>0);
fZ=@(t,P,Z) gamma*Rm*Z.*P.^2./(alpha^2+P.^2)-mu*Z;

u0 = 0.25; % for the T and B model this is in days
mydiff = 1e-3; 

%% Utilities for the numerical method
%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks = [0:N/2-1 0 -N/2+1:-1]*nyquist_freq;
ks2 = ks.*ks; ks3 = ks2.*ks; Dx = sqrt(-1)*ks;

% time stepping parameters
dt = 0.0005; numsteps = 200; numouts = 600
difffact = 1./(1+ks2*mydiff*dt);
Pf=fft(P);
Zf=fft(Z);
if LINE_PLOT==1
    figure(1)
    clf
    betterplots
    subplot(2,1,1)
    plot(x,P,'Color',rand(3,1,1))
    grid on
    ylabel('P')
    xlabel('x')
    subplot(2,1,1)
    plot(x,Z,'Color',rand(3,1,1))
    grid on
    ylabel('Z')
    xlabel('x')
    drawnow
end
t = 0;
Ps = zeros(numouts+1,N);
Zs = Ps;
Ps(1,:) = P;
Zs(1,:) = Z;
ts(1) = t;
%% Main time-stepping double loop
for ii = 1:numouts
    for jj = 1:numsteps
        % take the fft
        Pf = fft(P);
        Zf = fft(Z);
        % step forward in Fourier space
        reactPf = fft(fP(t,P,Z))*dt;
        reactZf = fft(fZ(t,P,Z))*dt;
        Pf = difffact.*(Pf-u0*dt*Dx.*Pf+reactPf);
        Zf = difffact.*(Zf-u0*dt*Dx.*Zf+reactZf);
        %transform back
        P = real(ifft(Pf));
        Z = real(ifft(Zf));      
        t = t+dt;
    end
    % make plots and/or save data
    if LINE_PLOT==1
        subplot(2,1,1)
        plot(x,P,'Color','k')
        axis([xmin xmax 0 100])
        grid on
        ylabel('P')
        title(['t = ' num2str(t,3)])
        subplot(2,1,2)
        plot(x,Z,'Color','r')
        axis([xmin xmax 0 25])
        grid on
        ylabel('Z')
        xlabel('x')
        drawnow
    end
    if SAVE_FIELDS==1
        Ps(ii+1,:) = P;
        Zs(ii+1,:) = Z;
        ts(ii+1) = t;
    end
end
if SAVE_FIELDS
    figure
    clf
    betterplots
    subplot(2,1,1)
    pcolor(x,ts,Ps),shading flat,colormap darkjet,colorbar
    ylabel('t')
    subplot(2,1,2)
    pcolor(x,ts,Zs),shading flat,colormap darkjet,colorbar
    xlabel('x')
    ylabel('t')
    subplot(2,1,1)
    text(-9,53,'(a)','Color','w')
    subplot(2,1,2)
    text(-9,53,'(b)','Color','w')
    
    figure
    clf'
    betterplots
    subplot(2,1,1)
    plot(x,Ps(101,:),'k',x,Ps(201,:),'r',x,Ps(301,:),'g',x,Ps(401,:),'m',x,Ps(501,:),'b')
    ylabel('P')
    legend('t=10','t=20','t=30','t=40','t=50')
    subplot(2,1,2)
    plot(x,Zs(101,:),'k',x,Zs(201,:),'r',x,Zs(301,:),'g',x,Zs(401,:),'m',x,Zs(501,:),'b')
    ylabel('Z')
    xlabel('x')
end