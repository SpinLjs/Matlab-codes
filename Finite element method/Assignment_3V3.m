%clear all;

% T=zeros(1,200);
% f_T=zeros(1,200);
% 
% for i2=1:length(T)
% 
% f=1e9+9e9/length(T)*(i2-1);   %Hz
f=40e9;
% f_T(i2)=f;

u0=4*pi*1e-7;   %N/A^2
e0=8.854187817e-12; %F/m  \

w=2*pi*f;
k0=w*sqrt(u0*e0);
a1=exp(j*0);


% z=[0 20e-3 21e-3 30e-3 31e-3 40e-3];
% p_z=[1 4 1.15 4 1; 1 1 1 1 1];
z=[0 20e-3 30e-3 50e-3];
p_z=[1 3-0.1j 1; 1 3-0.1j 1 ];

seg_N=400;
[z_seg, p_z_seg]=ass3_seg(z,p_z,seg_N);
er=p_z_seg(1,1:end);
ur=p_z_seg(2,1:end);
h_e=diff(z_seg);
z_e1=z_seg(1:end-1);
z_e2=z_seg(2:end);

Weight=[0.129484966168870 0.279705391489277 0.381830050505119 0.417959183673469 0.381830050505119 0.279705391489277 0.129484966168870];
Node=[-0.949107912342759 -0.741531185599394 -0.405845151377397 0 0.405845151377397 0.741531185599394 0.949107912342759];

A=kron(-1./ur./h_e,[1 -1;-1 1]);
B=kron(er.*h_e./6,[2 1;1 2]);

K_mid=A+k0^2*B;

K_diag0=[K_mid(1:4:end) 0]+[0 K_mid(4:4:end)];
K_diag1=K_mid(3:4:end);
K_diagM1=K_mid(2:4:end);
K=diag(K_diag0)+diag(K_diag1,1)+diag(K_diagM1,-1);

% 1.24
%-------------------try 1-------------------------------------------
coef_g=a1.*(1./ur-er)*k0^2;
% 


z_gq=Node.'*h_e/2+repmat((z_e1+z_e2)/2,length(Node),1);
g_temp1=(repmat(z_e2,length(Node),1)-z_gq).*exp(-j*k0*z_gq); 
g_temp2=(z_gq-repmat(z_e1,length(Node),1)).*exp(-j*k0*z_gq);
g_e1=Weight*g_temp1.*coef_g/2;
g_e2=Weight*g_temp2.*coef_g/2;

G=[g_e1 0]+[0 g_e2];
 
P=zeros(size(z_seg));
% H_inc_p=a1/sqrt(u0/e0)*exp(-j*k0*z_seg);
% P(2:end-1)=j*w*u0*(ur(1:end-1)-ur(2:end)).*H_inc_p(2:end-1);
% P(801)=   5.2258e+02 + 1.1773e+02i;
% P(1201)=1.6230e+02 + 5.1050e+02i;

%-----ABC, source at z->-infinity--------------------------
% 
K(1,1)=K(1,1)-j*k0/ur(1);                          % ABC
K(end,end)=K(end,end)-j*k0/ur(end);

%------------------------------------------------------------

AH=kron(-1./ur./h_e,[1 -1;-1 1]);
BH=kron(er.*h_e./6,[2 1;1 2]);
KH_mid=A+k0^2*B;
KH_diag0=[KH_mid(1:4:end) 0]+[0 KH_mid(4:4:end)];
KH_diag1=KH_mid(3:4:end);
KH_diagM1=KH_mid(2:4:end);
KH=diag(KH_diag0)+diag(KH_diag1,1)+diag(KH_diagM1,-1);
KH(1,1)=KH(1,1)-j*k0/er(1);                          % ABC
KH(end,end)=KH(end,end)-j*k0/er(end);

CE=diag([0 -j*w*u0*diff(ur) 0]);
CH=diag([0 -j*w*e0*diff(er) 0]);


%-----------------------------------------------
coef_gH=a1.*(1./er-ur)*k0^2;
% 


z_gqH=Node.'*h_e/2+repmat((z_e1+z_e2)/2,length(Node),1);
% z_gqH=Node.'*h_e+repmat(z_e1,length(Node),1);
g_temp1H=(repmat(z_e2,length(Node),1)-z_gqH).*exp(-j*k0*z_gqH); 
g_temp2H=(z_gqH-repmat(z_e1,length(Node),1)).*exp(-j*k0*z_gqH);
g_e1H=Weight*g_temp1H.*coef_gH/2;
g_e2H=Weight*g_temp2H.*coef_gH/2;
GH=[g_e1 0]+[0 g_e2];

E_inc_m=a1*exp(-j*k0*z_seg);
H_inc_m=a1/sqrt(u0/e0)*exp(-j*k0*z_seg);
RE=CE*H_inc_m.';
RH=CH*E_inc_m.';
L=K-CE*inv(KH)*CH;
R=G.'+RE+CE*inv(KH)*(GH.'+RH);
E_scat=inv(L)*R;



%-----------------------------------------------------
% F_scat=inv(K)*(G+P).';
% F_inc=a1*exp(-j*k0*z_seg);
% F=F_scat.' +F_inc;
% 
% figure(1)
% plot(z_seg,abs(F_scat));
% figure(2)
% plot(z_seg,abs(F));
% figure(3)
% plot(z_seg,phase(F_scat));
% figure(5);
% plot(z_seg,abs(G+P));
% figure(6);
% plot(z_seg,phase(G+P));
% figure(9);
% plot(abs([g_e1 0]),'b');
% hold on;
% plot(abs([0 g_e2]),'g');

figure(1)
plot(z_seg,abs(E_scat));
% figure(2)
% plot(z_seg,abs(E));
figure(3)
plot(z_seg,phase(E_scat));
% figure(5);
% plot(z_seg,abs(G+P));
% figure(6);
% plot(z_seg,phase(G+P));
% figure(9);
% plot(abs([g_e1 0]),'b');
% hold on;
% plot(abs([0 g_e2]),'g');
