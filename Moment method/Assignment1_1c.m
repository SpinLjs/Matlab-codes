clc;
clear all;
itv = 1e-3;             % Z axis point-to-point distance(all distance unit is m)  
step=5e-5;              % b-a, the width of quadrature strip
z=1e-3:itv:6e-1;        % To make the code be simple, I calculate Z on (0,0.6]
ZJ1=0.2;
ZJ2=0.32;
d=ZJ2-ZJ1;
Zc=(ZJ1+ZJ2)/2
ZM1=0.4;
ZM2=0.48;
f=5e9;                  % frequence
w=2*pi*f;               % angular frequence
u0=4*pi*1e-7;           % Vacuum permeability 
e0=8.854e-12;           % Vacuum permittivity F/m
k0=w*sqrt(e0*u0);
weight=[0.347854845 0.652145155 0.652145155 0.347854845];   % weight of 4-points
xi=[-0.861136312 -0.339981044 0.339981044 0.861136312];     % Xi values

Ex1=zeros(1,length(z));
Ex2=zeros(1,length(z));
Ex=zeros(1,length(z));
for i1=1:length(z);
    I1=0;               % For every i1, I1 equal to one strip's integration
    I2=0;               % I1 for electric current density, I2 for magnetic
    
    a=ZJ1;              % Ex1 from electric current density
    b=a+step;
    while(b<ZJ2);       % caluculate integration strip-by-strip
        x=(b-a)/2*xi+(b+a)/2; % express Z' with Xi
        f1=(b-a)/2*exp(-j*k0*abs(z(i1)-x)).*(cos(pi/d*(x-Zc))+0.1i*sin(pi/d*(x-Zc))); % g(x)to f(Xi);
        I1=I1+weight*f1.';   % product of weight vector and f(Xi(n)) vector;
        a=a+step;       % go to next integration strip 
        b=b+step;
    end
    Ex1(i1)=-w*u0/2/k0*I1; % get E field of electric current density;
    
    a=ZM1;                 % Ex2 from magnetc current density
    b=a+step;
    while(b<ZM2);          % caluculate integration strip-by-strip
        x=(b-a)/2*xi+(b+a)/2;   % express Z' with Xi
        if b<=z(i1);       % Gem's sign is decided by sign of (z-z')
        f2=(b-a)/2*exp(-j*k0*abs(z(i1)-x));     % g(x)to f(Xi);
        else
        f2=-(b-a)/2*exp(-j*k0*abs(z(i1)-x)) ; 
        end
        I2=I2+weight*f2';  % product of weight vector and f(Xi(n)) vector; 
        a=a+step;          % go to next integration strip
        b=b+step;
    end
    Ex2(i1)=sqrt(u0/e0)/2*I2;
    
end
Ex=Ex1+Ex2;             % total E field equal Ex1+Ex2

figure(1);                      % plot six graphs
title('Magnitude|Ex|');
plot(z(10:end),abs(Ex(10:end)));
xlabel('z/mm');
ylabel('Magnitude|Ex|');
figure(2);
title('Phase of Ex');
plot(z(10:end),phase(Ex(10:end)));
xlabel('z/mm');
ylabel('Phase of Ex');
figure(3);
title('Magnitude|Ex1|');
plot(z(10:end),abs(Ex1(10:end)));
xlabel('z/mm');
ylabel('Magnitude|Ex1|');
figure(4);
title('Phase of Ex1');
plot(z(10:end),phase(Ex1(10:end)));
xlabel('z/mm');
ylabel('Phase of Ex1');
figure(5);
title('Magnitude|Ex2|');
plot(z(10:end),abs(Ex2(10:end)));
xlabel('z/mm');
ylabel('Magnitude|Ex2|');figure(4);
figure(6)
title('Phase of Ex2');
plot(z(10:end),phase(Ex2(10:end)));
xlabel('z/mm');
ylabel('Phase of Ex2');


