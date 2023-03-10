%%FFT based method for Fisher's equation.
% Uses implicit scheme and compares to either an explicit scheme
% or a weaker diffusion
clear all; close all;
% plotting flags
LINE_PLOT = 1;
SAVE_FIELDS = 1;
%% set paramaters for the model here
nu = 1e-3; 
r = 5;
%nu=2e-1; % explicit goes unstable for this value of diffusivity

%% Utilities for the numerical method
%define a spatial grid
xmin = -10;
xmax = 10;
N = 2^12;
x = linspace(xmin,xmax,N);
dx = x(2)-x(1);

%make initial condition
u1 = 0.5*(1+tanh((x+5)/0.1))-0.5*(1+tanh((x-5)/0.1));
u2 = u1;%.*sin(2*pi*x/2);
u0 = u1;


%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks = [0:N/2-1 0 -N/2+1:-1]*nyquist_freq;
ks2 = ks.*ks; ks3 = ks2.*ks; ks4 = ks2.^2;

% define the filter
% Here I use hypervisocsity so the order should be even
% The parameter can be somewhat automated but here I tune by hand
filtord = 8;
hypervisc = 1e-14;
myfilt = 1./(1+hypervisc*ks.^filtord);

%time step and number of steps
t = 0;
dt = 4e-4;
impop = 1./(1+dt*nu*ks2);
impop2 = 1./(1+dt*0.1*nu*ks2);
%impop = 1./(1+dt*nu*ks4); % hyperviscosity

expop = -dt*ks2*nu;
numstps = 100;
numouts = 100;
if LINE_PLOT==1
    figure(1)
    clf
    % this sets the thick line plotting parameters useful for hard copies and journals
     set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
            'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
              'DefaultAxesFontWeight','bold');
    subplot(2,1,1)
    plot(x,u0,'k-')
    subplot(2,1,2)
    u0f = fft(u0);spec0=log10(u0f.*conj(u0f));
    plot(fftshift(ks),fftshift(spec0),'k.-')
end
us1=zeros(N,numouts+1);
us2=zeros(N,numouts+1);
us1(:,1) = u1;
us2(:,1) = u2;
ts(1) = 0;

%start with Euler time stepping
for ii = 1:numouts
   for jj = 1:numstps 
    t = t+dt;
    % base case
    u1 = real(ifft(impop.*(fft(u1 +dt*r*u1.*(1-u1)))));
    % manipulation (weaker diffusion)
    u2 = real(ifft(impop2.*(fft(u2 +dt*r*u2.*(1-u2)))));
    % manipulation (explicit time stepping)
    %u2 = real(ifft((fft(u2) +dt*r*u2.*(1-u2)+expop.*fft(u2))));
   end 
   if LINE_PLOT==1
       figure(1)
       clf
       subplot(2,1,1)
       plot(x,u0,'k-',x,u1,'b-',x,u2,'r-')
       grid on
       xlabel('x');
       ylabel('u');
       title(['time = ' num2str(t,2)]);
       axis([xmin xmax -0.25*max(max(u0)) 1.4*max(max(u0))])
       legend('initial','base case','manipulation','Location','South')
       subplot(2,1,2)
       u1f=fft(u1);spec1=log10(u1f.*conj(u1f));u2f=fft(u2);spec2=log10(u2f.*conj(u2f));
       plot(fftshift(ks),fftshift(spec0),'k.-',fftshift(ks),fftshift(spec1),'b.-',fftshift(ks),fftshift(spec2),'r.-')
       xlabel('k');
       ylabel('spectrum');
       grid on
       axis([-max(abs(ks(:))) max(abs(ks(:))) -35 8])
       drawnow 
   end
   if SAVE_FIELDS
       us1(:,ii+1) = u1;
       us2(:,ii+1) = u2;
       ts(ii+1) = t;
   end
end
if SAVE_FIELDS==1
    figure
    clf
    betterplots
    colormap darkjet
    subplot(2,1,1)
    pcolor(x,(0:numouts)*dt*numstps,us1'),shading flat
    ylabel('time')
    axis([3 7 0 max(ts(:))])
    subplot(2,1,2)
    pcolor(x,(0:numouts)*dt*numstps,us2'),shading flat
    ylabel('time')
    xlabel('x')
    axis([3 7 0 max(ts(:))])
end