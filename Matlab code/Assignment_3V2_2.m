%clear all;

T=zeros(1,200);
f_T=zeros(1,200);
S12_T=zeros(1,200);
P_T=zeros(1,200);
% 
for i2=1:length(T)
% 
f=1e9+39e9/length(T)*(i2-1);   %Hz
% f=10e9;


f_T(i2)=f;

u0=4*pi*1e-7;   %N/A^2
e0=8.854187817e-12; %F/m  \

w=2*pi*f;
k0=w*sqrt(u0*e0);
a1=exp(j*0);

% 
% z=[0 20e-3 21e-3 30e-3 31e-3 40e-3];
% p_z=[1 4 1.15 4 1; 1 1 1 1 1];
% z=[0 20e-3 28e-3 30e-3 38e-3 50e-3];
% p_z=[1 2-0.1j 40 1.22-1j 1; 1 2-0.1j 40 1.22-1j 1 ];
z=[0 20e-3 30e-3 50e-3];
p_z=[1 1-0.1j 1; 1 1-0.1j 1 ];

seg_N=60;
[z_seg, p_z_seg]=ass3_seg(z,p_z,seg_N);
er=p_z_seg(1,1:end);
ur=p_z_seg(2,1:end);
h_e=diff(z_seg);
z_e1=z_seg(1:end-1);
z_e2=z_seg(2:end);

Weight=[0.129484966168870 0.279705391489277 0.381830050505119 0.417959183673469 0.381830050505119 0.279705391489277 0.129484966168870];
Node=[-0.949107912342759 -0.741531185599394 -0.405845151377397 0 0.405845151377397 0.741531185599394 0.949107912342759];

A=kron(-1./ur./h_e,[1 -1;-1 1]);　　　　% Kronecker product
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
% z_gq=Node.'*h_e+repmat(z_e1,length(Node),1);
g_temp1=(repmat(z_e2,length(Node),1)-z_gq).*exp(-j*k0*z_gq); 
g_temp2=(z_gq-repmat(z_e1,length(Node),1)).*exp(-j*k0*z_gq);
g_e1=Weight*g_temp1.*coef_g/2;  
g_e2=Weight*g_temp2.*coef_g/2;

G=[g_e1 0]+[0 g_e2];
 
P=zeros(size(z_seg));
E_inc=a1*exp(-j*k0*z_seg);
P(2:end-1)=(1./ur(1:end-1)-1./ur(2:end))*-j*k0*a1.*exp(-j*k0*z_seg(2:end-1));   % inner boundary condition for scattered field

%-----ABC, source at z->-infinity--------------------------
% 
K(1,1)=K(1,1)-j*k0/ur(1);                          % ABC
K(end,end)=K(end,end)-j*k0/ur(end);

F_scat=inv(K)*(G+P).';
F_inc=a1*exp(-j*k0*z_seg);
F=F_scat.' +F_inc;



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





%-------------------------post process----------------------
%%
%----------S paramter------------------

z_in=10e-3;
z_out=38e-3;
[USLS,N_z_in]=min(abs(z_seg-z_in));
[USLS,N_z_out]=min(abs(z_seg-z_out));
S_11=F_scat(N_z_in)/F_inc(N_z_in);
S_21=F(N_z_out)/F_inc(N_z_in);

T(i2)=S_11;
S12_T(i2)=S_21;

plot(f_T,mag2db(abs(T)));
% 
% %%
% %---------- Power dissipation----------
    P_e_loss=zeros(size(er));
for i=1:length(z_seg)-1
    B_temp=B(:,2*i-1:2*i);
    P_e_loss(i)=w*e0/2*(F(i:i+1)*imag(B_temp)*(F(i:i+1)'));
end;
    p_loss=sum(P_e_loss);
    P_T(i2)=p_loss;
    
end
figure(21)
plot(f_T,mag2db(abs(T)));
figure(22)
plot(f_T,mag2db(abs(S12_T)));
figure(23)
plot(f_T,abs(P_T));
%     
%     
% %%    
% %----------- equivalent volume ------------
% 
% j_eq=j*w*e0.*([er 1]-1).*F;
% 
% 
% 
% %z_gq=Node.'*h_e/2+repmat((z_e1+z_e2)/2,length(Node),1);
% %z_gq=Node.'*h_e+repmat(z_e1,length(Node),1);
% z_ob=0:5e-6:40e-3;
% E_sc_ob=zeros(size(z_ob));
% M_sc_ob=zeros(size(z_ob));
% 
% E_inc=a1*exp(-j*k0*z_ob);
% H_inc=a1/sqrt(u0/e0)*exp(-j*k0*z_ob);
% 
% for i=1:length(z_ob) 
%     Psi_e_t1=(repmat(z_e2,length(Node),1)-z_gq).*exp(-j*k0*abs(z_ob(i)-z_gq)); 
%     Psi_e_t2=(z_gq-repmat(z_e1,length(Node),1)).*exp(-j*k0*abs(z_ob(i)-z_gq));
%     Psi_m_t1=(repmat(z_e2,length(Node),1)-z_gq).*(-1*sign(z_ob(i)-z_gq).*exp(-j*k0*abs(z_ob(i)-z_gq))); 
%     Psi_m_t2=(z_gq-repmat(z_e1,length(Node),1)).*(-1*sign(z_ob(i)-z_gq).*exp(-j*k0*abs(z_ob(i)-z_gq)));
%     Psi_e1=j/2/2/k0.*(Weight*Psi_e_t1);
%     Psi_e2=j/2/2/k0.*(Weight*Psi_e_t2);
%     Psi_m1=1/2/2.*(Weight*Psi_m_t1);
%     Psi_m2=1/2/2.*(Weight*Psi_m_t2);
%     E_seg=(er-1).*(F(1:end-1).*Psi_e1+F(2:end).*Psi_e2);
%     H_seg=(er-1).*(F(1:end-1).*Psi_m1+F(2:end).*Psi_m2);
%     E_sc_ob(i)=j*w*e0*sum(E_seg)*j*w*u0;
%     H_sc_ob(i)=j*w*e0*sum(H_seg);
%     
% end
% 
% 
% E_ob=E_inc+E_sc_ob;
% H_ob=H_inc+H_sc_ob;


% figure(4)
% plot(z_ob,abs(E_sc_ob));
% figure(5)
% plot(z_ob,abs(E_ob));
% figure(6)
% plot(z_ob,abs(H_sc_ob));
% figure(7)
% plot(z_ob,abs(H_ob));
%}
