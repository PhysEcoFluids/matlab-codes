%% This script solves the advection reaction diffusion equation in one dimension
% An FFT based method is used
%% Definitions and parameters
%define a spatial grid
xmin = -1e1;
xmax = 1e1;
N = 2^9;
x = linspace(xmin,xmax,N+1); x=x(1:end-1);
dx = x(2)-x(1);
% set flags here for different options for plotting
LINE_PLOT = 1;
SAVE_FIELD = 1;
%make initial condition
C = sech(((x+5)/1).^2)+sech(((x-5)/1).^2).*cos(2*pi*x/0.4);
C0 = C;
% advection and diffusion parameters
u0 = 1;
mydiff = 0.1;

%% Utilities for the numerical method
%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks = [0:N/2-1 0 -N/2+1:-1]*nyquist_freq;
ks2 = ks.*ks; ks3 = ks2.*ks; Dx = sqrt(-1)*ks;

% time stepping paramaters
dt=0.001;numsteps = 10;numouts = 100
% define the factor for implicit diffusion
difffact = 1./(1+ks2*mydiff*dt);
Cf = fft(C);
% if you want line plots make one of those for early times
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
%save initial conditions
t = 0;
Cs(1,:) = C;
ts(1) = t;
%% Main time stepping double loop
for ii = 1:numouts
    for jj = 1:numsteps
        Cf = difffact.*(Cf-u0*dt*Dx.*Cf+Cf*dt);
        t = t+dt;
    end
    % plot and/or save fields if flags are set for that
    C = real(ifft(Cf));
    if LINE_PLOT==1
        plot(x,C,'Color',rand(3,1,1))
        drawnow
    end
    if SAVE_FIELD==1
        Cs(ii+1,:)=C;
        ts(ii+1)=t;
    end
end
% Make a Hovmoller plots at the end
figure
clf
betterplots
pcolor(x,ts,Cs),shading flat,colormap darkjet,colorbar
xlabel('x')
ylabel('t')