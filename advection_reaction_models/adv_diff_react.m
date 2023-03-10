% number of grid points (total is N+1) and the grid
clear all,close all
N=64*4;
z=linspace(-1,1,N+1);
z=z(1:end-1);
Lx=1;Ly=1;
%[x,y]=meshgrid(2*z*L,z*L);
[x,y]=meshgrid(z*Lx,z*Ly);
% define three arrays for u now u in the past and u in the future
cn=exp(-20*(x.*x+y.*y));
kappa0=1e-3


dk=2*pi/(2*Lx); dl=2*pi/(2*Ly);
ksvec(1)=0; ksvec(N/2+1)=0;
for ii=2:(N/2)
   ksvec(ii)=ii-1;
   ksvec(N/2+ii)=-N/2 + ii -1;
end
lsvec=ksvec*dl;ksvec=ksvec*dk; 
[k,l]=meshgrid(ksvec,lsvec);
k2=k.*k; l2=l.*l; lap=-k2-l2;

 hypervisc=1e-16; filtord=8;
 myfilta=1./(1+hypervisc*(k.^filtord+l.^filtord));

psi=sin(2*pi*x/Lx).*cos(2*pi*y/Ly);
v=-real(ifft2(sqrt(-1)*k.*fft2(psi)));
u=real(ifft2(sqrt(-1)*l.*fft2(psi)));

dt=5e-5;
dt2=dt*dt;
numsteps=2000;
numouts=8;

imdiff=1./(1-dt*kappa0*lap);
imdiff2=1./(1-10*dt*kappa0*lap);
myfilt=ones(size(imdiff));
% kill the Nyquist filter
myfilt(N/2-2:N/2+4,:)=0;myfilt(:,N/2-2:N/2+4)=0;
myfilt(:,N/2-2:N/2+4)=0;myfilt(N/2-2:N/2+4,:)=0;


cn2=cn;
figure(1)
clf

pcolor(x,y,cn); shading interp
xlabel('x','fontweight','bold','fontsize',12);
ylabel('y','fontweight','bold','fontsize',12);
%caxis([-0.2,0.2])
axis([-Lx Lx -Ly Ly])
colorbar
drawnow
t=0
figure
clf
colormap jet
%subplot(3,3,1)
pcolor(x,y,cn); shading interp
for ii=1:numouts
 for jj=1:numsteps
   % timestep using fft2 euler for now
   t=t+dt;
   dummy=fft2(cn);
   dummy2=fft2(cn2);
   cnx=real(ifft2(sqrt(-1)*k.*dummy));
   cny=real(ifft2(sqrt(-1)*l.*dummy));
   cn2x=real(ifft2(sqrt(-1)*k.*dummy2));
   cn2y=real(ifft2(sqrt(-1)*l.*dummy2));
   clap=real(ifft2(lap.*dummy));
   source=cn.*(5-cn);
   source2=cn2.*(5-cn2);
   cn=real(ifft2(myfilt.*imdiff.*fft2(cn+dt*(-u.*cnx-v.*cny+source))));
   cn2=real(ifft2(myfilt.*imdiff2.*fft2(cn2+dt*(-u.*cn2x-v.*cn2y+source2))));
   %% implicit diffusion

 end
 ii+1
 %subplot(3,3,ii+1)
 figure
colormap jet
subplot(2,2,1)
 pcolor(x,y,cn); shading interp
 title(num2str(t,3))
 axis([-Lx Lx -Ly Ly])
 subplot(2,2,2)
 pcolor(x,y,cn2); shading interp
 title(num2str(t,3))
 axis([-Lx Lx -Ly Ly])
 drawnow


end

