%---------------------------(a)--------------------------------------------
% x_min=-10;
% d_x=0.5;
% x_max=10;
% y_min=-10;
% d_y=0.5;
% y_max=10;
% 
% [x,y]=meshgrid(x_min:d_x:x_max,y_min:d_y:y_max);
% Px=zeros(size(x));
% Py=Px;
% Px=x.*y;
% Py=y;
% 
% h=max([max(abs(Px)) max(abs(Py))])./min([d_x d_y]);
% 
% my_quiver2(x,y,Px,Py);

%---------------------------(b)--------------------------------------------
% clear all;
% f=9.5e9;   %Hz
% a=22.86e-3;    %m
% t=0;
% z=0;
% 
% u0=4*pi*1e-7;   %N/A^2
% e0=8.854187817e-12; %F/m  
% w=2*pi*f;
% k=w*sqrt(u0*e0);
% m=1;
% k_c=sqrt((m*pi/a).^2);
% beta=sqrt(k.^2-k_c.^2);
% 
% x_min=0;
% d_x=a./100;
% x_max=a;
% y_min=0;
% y_max=10.16e-3;
% d_y=y_max./100;
% 
% [x,y]=meshgrid(x_min:d_x:x_max,y_min:d_y:y_max);
% Px=zeros(size(x));
% Py=Px;
% 
% Py=real(sin(pi*x/a)*exp(-j*beta*z)*exp(j*w*t));
% 
% h=max([max(abs(Px)) max(abs(Py))])./min([d_x d_y]);
% 
% %quiver(x,y,Px,Py);
% my_quiver2(x,y,Px,Py,h);

%---------------------------(c)--------------------------------------------
clear all;
f=9.5e9;   %Hz
y=5.08e-3;
a=22.86e-3;

u0=4*pi*1e-7;   %N/A^2
e0=8.854187817e-12; %F/m  
w=2*pi*f;
k=w*sqrt(u0*e0);
m=1;
k_c=sqrt((m*pi/a).^2);
beta=sqrt(k.^2-k_c.^2);

x_min=0;
x_max=22.86e-3;
d_x=(x_max-x_min)./800;
z_min=0;
z_max=100e-3;
d_z=(z_max-z_min)./1000;

[x,z]=meshgrid(x_min:d_x:x_max,z_min:d_z:z_max);
Px=zeros(size(x));
Py=Px;

t=0;
Py=real(sin(pi*x/a).*exp(-j*beta*z).*exp(j*w*t));

S1=mesh(z(:,1),x(1,:),Py.');
view(2);
colorbar;
caxis([-1,1]);

for t=5e-12:5e-12:1e-10;
    Py=real(sin(pi*x/a).*exp(-j*beta*z).*exp(j*w*t));
    set(S1,'ZData',Py.')
    drawnow;
    pause(0.05);
end






