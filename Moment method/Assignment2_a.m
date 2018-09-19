clear all;
clc;
zitv = 1e-3;
deltaz=zitv;            % 
zin=0;
f=1e9;                  % frequence
loss_tan=0.02;
er=[1,4.17*(1-j*loss_tan),1];       % epsilonr values for three slab
% er=[4 1.15 4];
d=[0 6e-1 0];     % width for three slabs
dlu=[1e-1 1e-1];        % distances between slab and observe points
E0=1;                   % 

zout=sum(d)+sum(dlu);             
z=zin:zitv:zout;        % 
epsilonr=[er(1)*ones(1,int16(d(1)/deltaz)) er(2)*ones(1,int16(d(2)/deltaz)) er(3)*ones(1,int16(d(3)/deltaz))];

w=2*pi*f;               % angular frequence
u0=4*pi*1e-7;           % Vacuum permeability 
e0=8.854e-12;           % Vacuum permittivity F/m
k0=w*sqrt(e0*u0);       % k0
eta0=sqrt(u0/e0);       % eta0
% Gauss-Legendre Quadrature(4 points)
weight=[0.347854845 0.652145155 0.652145155 0.347854845];   % weight of 4-points
xi=[-0.861136312 -0.339981044 0.339981044 0.861136312];     % Xi values

zc_1=z(int16(dlu(1)/zitv+1):int16(dlu(1)/zitv)+int16(sum(d)/zitv));
idx=find(epsilonr~=1);
zc=zc_1(idx)+deltaz/2;
epslr_zc=epsilonr(idx);

I_total=zeros(1,length(zc_1));
Zmn=zeros(length(zc),length(zc));
V=zeros(length(zc),1);
I=zeros(length(zc),1);
for(m=1:length(zc))
    for(n=1:length(zc))
        if(zc(m)==zc(n))
            Zmn(m,n)=1/(j*w*e0*(epslr_zc(m)-1))-j*eta0/k0*(1-exp(-j*k0*deltaz/2));
        else
            Zmn(m,n)=eta0/k0*sin(k0*deltaz/2)*exp(-j*k0*abs(zc(m)-zc(n)));
        end
    end
    V(m)=E0*exp(-j*k0*zc(m));
end
I=inv(Zmn)*V;
I_total(idx)=I;


Ex_scat=zeros(1,length(z));    % pre-difine the scale of E fields
Ex=zeros(1,length(z));
for i1=1:length(z);
    Itg1=0;               % For every i1, I1 equal to one strip's integration
    for i2=1:length(I_total)       % caluculate integration strip-by-strip
        if(I_total(i2)~=0)
            z1n=zc_1(i2);              % Ex1 from electric current density
            z2n=z1n+deltaz;
            x=(z2n-z1n)/2*xi+(z2n+z1n)/2; % express Z' with Xi
            f1=(z2n-z1n)/2*exp(-j*k0*abs(z(i1)-x))*I_total(i2); % g(x)to f(Xi);
            Itg1=Itg1+weight*f1.';   % product of weight vector and f(Xi(n)) vector;
        end
    end
    Ex_scat(i1)=-eta0/2*Itg1; % get E field of electric current density;
end

Ex_inc=E0*exp(-j*k0.*z);
Ex=Ex_scat+Ex_inc;

Reflect_coe=abs(Ex_scat(1)/Ex_inc(1))
Trans_coe=abs(Ex(end)/Ex_inc(1))

Ex_zc=Ex(ceil(zc./zitv));
I1=0;
for i1=1:length(zc)
    I1=I1+w*-imag(e0*epslr_zc(i1))*abs(Ex_zc(i1)).^2*deltaz;
end
P_loss=I1
P_in=1/eta0*abs(Ex_inc(1)).^2;
P_out=1/eta0*abs(Ex(end)).^2+1/eta0*abs(Ex_scat(1)).^2;
P_loss=P_in-P_out

figure(1);                      % plot six graphs
title('Magnitude|Ex|');
plot(z,abs(Ex));
xlabel('z/m');
ylabel('Magnitude|Ex|');
figure(2);
title('Phase of Ex');
plot(z,angle(Ex));
xlabel('z/m');
ylabel('Phase of Ex');
figure(3);
title('Jx_eq');
plot(zc_1,abs(I_total));
xlabel('z/m');
ylabel('Jx_eq');
figure(4);                     
title('Magnitude|Ex_inc|');
plot(z,abs(Ex_inc));
xlabel('z/m');
ylabel('Magnitude|Ex_inc|');
figure(5);                     
title('Magnitude|Ex_scat|');
plot(z,abs(Ex_scat));
xlabel('z/m');
ylabel('Magnitude|Ex_scat|');


