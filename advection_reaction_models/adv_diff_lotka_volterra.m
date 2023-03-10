%% This script solves the Lotke Volterra model with advection and diffusion
clear all
close all
%% Definitions and parameters
%define a spatial grid
xmin = -1e1;
xmax = 1e1;
N = 2^9;
x = linspace(xmin,xmax,N+1); x=x(1:end-1);
dx = x(2)-x(1);
LINE_PLOT = 1;
SAVE_FIELDS = 1;
%make initial condition
C = 0.3+0.4*(sech(((x+5)/1).^2)+sech(((x-5)/1).^2).*cos(2*pi*x/0.4));
F = 0.05+0.1*rand(size(C));
C0 = C;
u0 = 0.25; % for the longer integration time 0.25 is a good choice for figures
%mydiff = 1e-2; % turn this down to 1e-3 to see some effect of packets
mydiff = 1e-3; % turn this up to 1e-2 to smooth out the effect of packets
mya = 0.2;

%% Utilities for the numerical method
%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks = [0:N/2-1 0 -N/2+1:-1]*nyquist_freq;
ks2 = ks.*ks; ks3 = ks2.*ks; Dx = sqrt(-1)*ks;

dt = 0.0005; numsteps = 200; numouts = 400
difffact = 1./(1+ks2*mydiff*dt);
Cf = fft(C);
Ff = fft(C);
if LINE_PLOT==1
    figure(1)
    clf
    betterplots
    plot(x,C,'Color',rand(3,1,1))
    hold on
    grid on
    ylabel('C')
    xlabel('x')
    drawnow
end
t = 0;
Cs = zeros(numouts+1,N);
Fs = Cs;
Cs(1,:) = C;
Fs(1,:) = F;
ts(1) = t;
%% Main time-stepping double loop
for ii = 1:numouts
    for jj = 1:numsteps
        % We can do most timestepping in Fourier space
        % except for the nonlinear term so we take care of that here
        nonlintermf = fft(real(ifft(Cf)).*real(ifft(Ff)));
        % reactions (this might be nicer as inline functions)
        reactCf = (mya*Cf-nonlintermf)*dt;
        reactFf = (-Ff+nonlintermf)*dt;
        % step forward
        Cf = difffact.*(Cf-u0*dt*Dx.*Cf+reactCf);
        Ff = difffact.*(Ff-u0*dt*Dx.*Ff+reactFf);
        t = t+dt;
    end
    % make plots and/or save data
    C = real(ifft(Cf));
    F = real(ifft(Ff));
    if LINE_PLOT==1
        subplot(2,1,1)
        plot(x,C,'Color','k')
        axis([xmin xmax -0.1 2.5])
        grid on
        ylabel('C')
        title(['t = ' num2str(t,3)])
        subplot(2,1,2)
        plot(x,F,'Color','r')
        axis([xmin xmax -0.1 1.1])
        grid on
        ylabel('F')
        xlabel('x')
        drawnow
    end
    if SAVE_FIELDS==1
        Cs(ii+1,:) = C;
        Fs(ii+1,:) = F;
        ts(ii+1) = t;
    end
end
if SAVE_FIELDS
    figure
    clf
    betterplots
    subplot(2,1,1)
    pcolor(x,ts,Cs),shading flat,colormap darkjet,colorbar
    ylabel('t')
    subplot(2,1,2)
    pcolor(x,ts,Fs),shading flat,colormap darkjet,colorbar
    xlabel('x')
    ylabel('t')
    subplot(2,1,1)
    text(-9,37,'(a)','Color','w')
    subplot(2,1,2)
    text(-9,37,'(b)','Color','w')
end