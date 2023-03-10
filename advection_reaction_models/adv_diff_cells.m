%% This script solves the 2D advection and diffusion problem with a multi-cell flow
clear all,close all
%% Definitions and parameters
N = 128;
z = linspace(-1,1,N+1);
z = z(1:end-1);
L = 1;
[x,y] = meshgrid(z*L,z);
% define three arrays for u now u in the past and u in the future
cn = exp(-((x-0.5*L).^8+(y+0.85).^8)/(0.1^8));cn0=cn;
kappa0 = 1e-2

%% Utilities for the numerical method
dk = 2*pi/(2*L); dl = 2*pi/(2*1);
ksvec(1) = 0; ksvec(N/2+1) = 0;
for ii=2:(N/2)
   ksvec(ii) = ii-1;
   ksvec(N/2+ii) = -N/2 + ii -1;
end
ksvec = ksvec*dk; lsvec = ksvec*dl;
[k,l] = meshgrid(ksvec,lsvec);
k2 = k.*k; l2 = l.*l; lap = -k2-l2;

% The background flow
psi = sin(pi*x/L).*cos(0.5*pi*y/1);
u = -0.5*pi*sin(pi*x/L).*sin(0.5*pi*y/1);
v = -(pi/L)*cos(pi*x/L).*cos(0.5*pi*y/1);


% time stepping parameters
dt = 1e-5;
numsteps = 4*3600;
numouts = 10;

% implicit diffusion in Fourier space
imdiff = 1./(1-0.1*dt*kappa0*lap);
imdiff2 = 1./(1-dt*kappa0*lap);

cn2 = cn;
figure(1)
clf
pcolor(x,y,cn); shading interp
xlabel('x','fontweight','bold','fontsize',12);
ylabel('y','fontweight','bold','fontsize',12);
%caxis([-0.2,0.2])
axis([-L L -1 1])
title('Initial state')
colorbar
drawnow
t=0
%% main time stepping double loop
for ii = 1:numouts
 for jj = 1:numsteps
   % timestep using fft2 euler for now
   t = t+dt;
   dummy = fft2(cn);
   dummy2 = fft2(cn2);
   cnx = real(ifft2(sqrt(-1)*k.*dummy));
   cny = real(ifft2(sqrt(-1)*l.*dummy));
   cn2x = real(ifft2(sqrt(-1)*k.*dummy2));
   cn2y = real(ifft2(sqrt(-1)*l.*dummy2));
   % These two lines would add a Fisher-type source
   source = 0*cn.*(5-cn);
   source2 = 0*cn2.*(5-cn2);
   cn = real(ifft2(imdiff.*fft2(cn+dt*(-u.*cnx-v.*cny+source))));
   cn2 = real(ifft2(imdiff2.*fft2(cn2+dt*(-u.*cn2x-v.*cn2y+source2))));
   %% implicit diffusion
 end
 
 figure
 clf
 colormap darkjet
 subplot(3,1,1)
 pcolor(x,y,psi); shading interp
 hold on
 contour(x,y,cn0,10,'k')
 hold off
 ylabel('y','fontweight','bold','fontsize',12);
 axis([-L L -1 1])
 title(['Streamfunction (top), low diffusivity (middle), higher diffusivity (bot) at t = ' num2str(t,3)])
 subplot(3,1,2)
 pcolor(x,y,cn); shading interp,caxis([0 1]*0.1)
 ylabel('y','fontweight','bold','fontsize',12);
 axis([-L L -1 1])
 subplot(3,1,3)
 pcolor(x,y,cn2); shading interp,caxis([0 1]*0.1)
 xlabel('x','fontweight','bold','fontsize',12);
 ylabel('y','fontweight','bold','fontsize',12);
 axis([-L L -1 1])
 drawnow
end