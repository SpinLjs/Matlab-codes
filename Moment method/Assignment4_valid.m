% Including customed function MyBesselh
clear all;
clc;

f=15e9;                 % frequence/ Hz
w=2*pi*f;               % angular frequence
u0=4*pi*1e-7;           % Vacuum permeability 
e0=8.854e-12;           % Vacuum permittivity F/m
k=w*sqrt(e0*u0);        % k0

L=6e-1;
b=12e-3; %1/2*1./sqrt(u0*e0)/f;
Za=0;
Zb=Za+L;
Ya=0;
Yb=Ya+b;
LS=[-5e-2 6e-3];
J0=1;

zitvob=1e-3;
yitvob=1e-3;
z=0e-1:zitvob:6e-1;  % observe points;
y=LS(2);
z_ob=repmat(z,length(y),1);
y_ob=repmat(y',1,length(z));

E_last=zeros(1,length(z_ob));

N=[1:69 70:5:145 150:10:690 700:100:2000];
Err=zeros(1,length(N));
for i_t=1:length(N);

zitv = L/N(i_t);
yitv = b/N(i_t);

% Gauss-MRW Quadrature(5 points)
nod=[0.56522282050801e-2 0.73403717426523e-1 0.284957404462558 0.619482264084778 0.915758083004698];    % 5 points G-MRW 
weight=[0.210469457918546e-1 0.130705540744447 0.289702301671314 0.350220370120399 0.208324841671986];  % weight

PEC_1=[Za:zitv:Zb-zitv];
PEC1=[PEC_1;Ya*ones(1,length(PEC_1));PEC_1+zitv;Ya*ones(1,length(PEC_1))];
PEC_2=[Za:zitv:Zb-zitv];
PEC2=[PEC_2;Yb*ones(1,length(PEC_2));PEC_2+zitv;Yb*ones(1,length(PEC_2))];
% PEC_3=[0:yitv:b-yitv];
% PEC3=[0*ones(1,length(PEC_3));PEC_3;0*ones(1,length(PEC_3));PEC_3+yitv];
% PEC_4=[0:yitv:b-yitv];
% PEC4=[Zb*ones(1,length(PEC_4));PEC_4;Zb*ones(1,length(PEC_4));PEC_4+yitv];
PEC=[PEC1 PEC2];
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

E_inc=-w*u0*J0/4*MyBesselh(z_ob,LS(1),LS(1),y_ob,LS(2),LS(2),k);

Az=zeros([length(y) length(z)]);
Azn_temp=zeros([length(y) length(z) length(nod)]);
for n=1:length(I)
    for i=1:length(nod)
       Azn_temp(:,:,i)=MyBesselh(z_ob,PEC(1,n),PEC(3,n),y_ob,PEC(2,n),PEC(4,n),k,nod(i));
       Azn_temp(:,:,i)=Azn_temp(:,:,i)*weight(i);
    end
    Az=Az+I(n)*sum(Azn_temp,3)*hn(n);
end
Az=Az*(u0/4j);
E_scat=-j*w*Az;
E_total=E_inc+E_scat;

Err(i_t)=sqrt(sum(abs(abs(E_total)-E_last).^2)/length(z));
E_last=abs(E_total);
end

step=ones(1,length(N));
for i=2:length(N)
step(i)=N(i)-N(i-1);
end
Err2=Err./step;

plot(N,10*log10(Err2));
xlabel('segments number per PEC plate N');
yabel('self-define parameter ERR/ dB');

