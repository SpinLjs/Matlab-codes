% Including customed function MyBesselh
clear all;
clc;

f=15e9;                 % frequence/ Hz
w=2*pi*f;               % angular frequence
u0=4*pi*1e-7;           % Vacuum permeability 
e0=8.854e-12;           % Vacuum permittivity F/m
k=w*sqrt(e0*u0);       % k0
n_TE=1;

r=3000;
Q=-pi:2*pi/360:pi; 
[r1,Q1]=meshgrid(r,Q);
z_ob=r1.*cos(Q1);
y_ob=r1.*sin(Q1);
z=zeros(1,length(r));
y=zeros(1,length(Q));



L=6e-1;
b=12e-3; %1/sqrt(u0*e0)/f/2;
Za=0;
Zb=Za+L;
Ya=0;
Yb=Ya+b;

k_c=n_TE*pi/b;      % cut-off wavenumbers of TEn
f_c=n_TE/(2*b*sqrt(u0*e0));     % cut-off frequencies
beta=sqrt(k^2-k_c^2);       % propagation constant         
lambda_g=2*pi/beta; % guide wavelength
d=lambda_g/4;

LS=[d b/2];
J0=1;

r_ls=sqrt(LS(1)^2+LS(2)^2);
Q_ls=atan(LS(2)/LS(1));

% zitvob=1e-3;
% yitvob=5e-4;
% z=-1e-1:zitvob:8e-1;  % observe points;
% y=-3e-2:yitvob:5e-2;
% z_ob=repmat(z,length(y),1);
% y_ob=repmat(y',1,length(z));

zitv = 5e-4;
yitv = 5e-5;

% Gauss-MRW Quadrature(5 points)
nod=[0.56522282050801e-2 0.73403717426523e-1 0.284957404462558 0.619482264084778 0.915758083004698];    % 5 points G-MRW 
weight=[0.210469457918546e-1 0.130705540744447 0.289702301671314 0.350220370120399 0.208324841671986];  % weight

PEC_1=[Za:zitv:Zb-zitv];
PEC1=[PEC_1;Ya*ones(1,length(PEC_1));PEC_1+zitv;Ya*ones(1,length(PEC_1))];
PEC_2=[Za:zitv:Zb-zitv];
PEC2=[PEC_2;Yb*ones(1,length(PEC_2));PEC_2+zitv;Yb*ones(1,length(PEC_2))];
PEC_3=[0:yitv:b-yitv];
PEC3=[0*ones(1,length(PEC_3));PEC_3;0*ones(1,length(PEC_3));PEC_3+yitv];
% PEC_4=[0:yitv:b-yitv];
% PEC4=[Zb*ones(1,length(PEC_4));PEC_4;Zb*ones(1,length(PEC_4));PEC_4+yitv];
PEC=[PEC1 PEC2 PEC3];
PEC_ct=[(PEC(1,:)+PEC(3,:))/2;(PEC(2,:)+PEC(4,:))/2];  %PEC center of segement

Vm=w*u0*J0/4*MyBesselh(PEC_ct(1,:),LS(1),LS(1),PEC_ct(2,:),LS(2),LS(2),k);

zm=repmat(PEC_ct(1,:).',1,length(PEC_ct(1,:)));
zn_1=repmat(PEC(1,:),length(PEC(1,:)),1);
zn_2=repmat(PEC(3,:),length(PEC(3,:)),1);
ym=repmat(PEC_ct(2,:).',1,length(PEC(2,:)));
yn_1=repmat(PEC(2,:),length(PEC(2,:)),1);
yn_2=repmat(PEC(4,:),length(PEC(4,:)),1);

Zmn=zeros(length(PEC_ct(1,:)));
Zmn_qrad=zeros([size(Zmn) length(nod)]);
hn=sqrt((PEC(3,:)-PEC(1,:)).^2+(PEC(4,:)-PEC(2,:)).^2);
hn_Z=repmat(hn,length(hn),1);
for i=1:length(nod)
   Zmn_qrad(:,:,i)=MyBesselh(zm,zn_1,zn_2,ym,yn_1,yn_2,k,nod(i));
   Zmn_qrad(:,:,i)=Zmn_qrad(:,:,i).*(1-diag(ones(1,length(Zmn))));
   xi0=1/2;
   Zmn_diag_1=xi0*MyBesselh(PEC_ct(1,:),PEC(1,:),PEC(3,:),PEC_ct(2,:),PEC(2,:),PEC(4,:),k,(1-nod(i))*xi0);
   Zmn_diag_2=(1-xi0)*MyBesselh(PEC_ct(1,:),PEC(1,:),PEC(3,:),PEC_ct(2,:),PEC(2,:),PEC(4,:),k,xi0+(1-xi0)*nod(i));
   Zmn_qrad(:,:,i)=Zmn_qrad(:,:,i)+diag(Zmn_diag_1+Zmn_diag_2);
   Zmn_qrad(:,:,i)=Zmn_qrad(:,:,i)*weight(i);
end
Zmn=sum(Zmn_qrad,3).*hn_Z*(-w*u0/4);

I=(inv(Zmn)*Vm.');

F_inc=-j*w*u0*J0/sqrt(8*j*pi*k).*exp(j*k*r_ls.*cos(Q1-Q_ls)); %*exp(-j*k*r1)./sqrt(r1)

Az=zeros([length(y) length(z)]);
Azn_temp=zeros([length(y) length(z) length(nod)]);
for n=1:length(I)
    for i=1:length(nod)
       Azn_temp(:,:,i)=exp(j*k*nod(i)*(PEC(3,n)-PEC(1,n)).*cos(Q1)+j*k*nod(i)*(PEC(4,n)-PEC(2,n)).*sin(Q1));
       Azn_temp(:,:,i)=Azn_temp(:,:,i)*weight(i);
    end
    Az=Az+I(n)*sum(Azn_temp,3)*hn(n).*exp(j*k*PEC(1,n)*cos(Q1)+j*k*PEC(2,n)*sin(Q1));
end
Az=Az*(u0/sqrt(8*j*pi*k));
F_scat=-j*w*Az;
F_total=F_inc+F_scat;


figure(1);
% mesh(z_ob,y_ob,10*log10(abs(F_total)));
mesh(z_ob,y_ob,abs(F_total));
figure(2);
mesh(z_ob,y_ob,abs(F_scat));
figure(3);
mesh(z_ob,y_ob,abs(F_inc));

figure(1);
% F_dB=10*log10(abs(F_total));
% F_dB(find(F_dB<0))=0;
F_dB=10*log10(abs(F_total)/max(abs(F_total)));
polar(Q.',F_dB);
figure(2);
F_dBscat=10*log10(abs(F_scat));
F_dBscat(find(F_dBscat<0))=0;
polar(Q.',F_dBscat);
figure(3);
F_dBinc=10*log10(abs(F_inc));
F_dBinc(find(F_dBinc<0))=0;
polar(Q.',F_dBinc);

