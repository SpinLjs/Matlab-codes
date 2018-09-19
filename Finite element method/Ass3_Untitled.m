clear all;

f=50e9;   %Hz
u0=4*pi*1e-7;   %N/A^2
e0=8.854187817e-12; %F/m  
w=2*pi*f;
k0=w*sqrt(u0*e0);
a1=exp(j*0);


z=[0 20e-3 21e-3 30e-3 31e-3 40e-3];
p_z=[1 1.15 4 1.15 1; 1 1 1 1 1];
seg_N=600;
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
coef_g=a1*k0^2.*(1./ur-er);
z_gq=Node.'*h_e/2+repmat((z_e1+z_e2)/2,length(Node),1);
z_gq=Node.'*h_e+repmat(z_e1,length(Node),1);
g_temp1=(repmat(z_e2,length(Node),1)-z_gq).*exp(-j*k0*z_gq); 
g_temp2=(z_gq-repmat(z_e1,length(Node),1)).*exp(-j*k0*z_gq);
g_e1=Weight*g_temp1;
g_e2=Weight*g_temp2;
g_e=repmat(coef_g,2,1).*[g_e1;g_e2].*1/2;

G=zeros(size(z_seg));
G=[g_e(1,:) 0]+[0 g_e(2,:)];

%-------------------try 2-------------------------------------------
% coef_g=a1*k0^2.*(1./ur-er);
% 
% Ng_seg=8;
% z_gseg=My_g_seg(z_seg,Ng_seg);
% h_gei=diff(z_gseg);
% z_gq=Node.'*h_gei+repmat(z_gseg(1:end-1),length(Node),1);
% % z_gq=Node.'*h_gei/2+repmat((z_gseg(1:end-1)+z_gseg(2:end))/2,length(Node),1);
% g_temp1=(repmat(z_gseg(2:end),length(Node),1)-z_gq).*exp(-j*k0*z_gq); 
% g_temp2=(z_gq-repmat(z_gseg(1:end-1),length(Node),1)).*exp(-j*k0*z_gq);
% g_e1=Weight*g_temp1;
% g_e2=Weight*g_temp2;
% g_eWE1=sum(reshape(g_e1,Ng_seg,length(g_e1)/Ng_seg));
% g_eWE2=sum(reshape(g_e2,Ng_seg,length(g_e1)/Ng_seg));
% g_e=repmat(coef_g,2,1).*[g_eWE1;g_eWE2].*1/2;
% 
% G=zeros(size(z_seg));
% G=[g_e(1,:) 0]+[0 g_e(2,:)];


P=zeros(size(z_seg));
%-----ABC, source at z->-infinity--------------------------
% P(1)=-j*2*k0*a1*exp(-j*k0*z_seg(1))/er(1);        % 1.21

K(1,1)=K(1,1)+j*k0/ur(1);                          % ABC
K(end,end)=K(end,end)-j*k0/ur(end);


F=inv(K)*((G+P).');

F_inc=a1*exp(-j*k0*z_seg);
F_tot=F.'+F_inc;

figure(1)
plot(z_seg,abs(F));
figure(2)
plot(z_seg,abs(F_tot));

%-------------------------post process----------------------








