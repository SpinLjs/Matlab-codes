clear all;

u0=4*pi*1e-7;   %N/A^2
e0=8.854187817e-12; %F/m 
%----- seven points Gauss-legender quaduature -------
Weight=[0.129484966168870 0.279705391489277 0.381830050505119 0.417959183673469 0.381830050505119 0.279705391489277 0.129484966168870];
Node=[-0.949107912342759 -0.741531185599394 -0.405845151377397 0 0.405845151377397 0.741531185599394 0.949107912342759];

a1=exp(j*0);    % Magnitude of input E field
%
%% -------------------- frequency sweep -------------------------------
Swp=100;                          
f_swp=1e9:9e9/(Swp-1):10e9;
S11e_T=zeros(1,Swp);
S21e_T=zeros(1,Swp);
S11m_T=zeros(1,Swp);
S21m_T=zeros(1,Swp);
P_loss_T=zeros(1,Swp);
% 
for i2=1:Swp
% 
f=f_swp(i2);
%}

% f=40e9;         % frequence
w=2*pi*f;
k0=w*sqrt(u0*e0);

%{
%% -------------------- convergency sweep -------------------------------
Swp=90;
N_swp =[1:49 50:5:50+(Swp-50)*5];                        
S11e_T=zeros(1,Swp);
S21e_T=zeros(1,Swp);
S11m_T=zeros(1,Swp);
S21m_T=zeros(1,Swp);
P_loss_T=zeros(1,Swp);
% 
for i2=1:Swp
%}



%% -------------- original setting for structure -------
z=[0 20e-3 21e-3 30e-3 31e-3 40e-3];
p_z=[1 1.15 4 1.15 1; 1 1 1 1 1 ];
% z=[0 15e-3 22e-3 30e-3 34e-3 37e-3 40e-3 50e-3];
% p_z=[1 2-0.1j 1 6 2.75-0.9j 1 1; 1 2-0.1j 1 1 3.2-0.8j 8 1 ];
% z=[0 20e-3 30e-3 50e-3];        
% p_z=[1 3-0.1j 1; 1 3-0.1j 1 ];  % first row: epsilon_r, second row: mu_r 

%% ----------- pre-processing ----------------------
seg_N=25;

[z_seg, p_z_seg]=ass3_seg(z,p_z,seg_N); % My costom mesh function ass3_seg()
er=p_z_seg(1,1:end);            % meshed epsilon_r
ur=p_z_seg(2,1:end);            % meshed mu_r
h_e=diff(z_seg);                % width of each element
z_e1=z_seg(1:end-1);
z_e2=z_seg(2:end);

Weight=[0.129484966168870 0.279705391489277 0.381830050505119 0.417959183673469 0.381830050505119 0.279705391489277 0.129484966168870];
Node=[-0.949107912342759 -0.741531185599394 -0.405845151377397 0 0.405845151377397 0.741531185599394 0.949107912342759];

%% ---------------------------- A B K -------------------------------
%------------------------------- for E ---------------------------
A=kron(-1./ur./h_e,[1 -1;-1 1]);    % Kronecker product
B=kron(er.*h_e./6,[2 1;1 2]);

K_mid=A+k0^2*B;

K_diag0=[K_mid(1:4:end) 0]+[0 K_mid(4:4:end)];
K_diag1=K_mid(3:4:end);
K_diagM1=K_mid(2:4:end);
K=diag(K_diag0)+diag(K_diag1,1)+diag(K_diagM1,-1);

%------------------------------------for H ------------------------

Am=kron(-1./er./h_e,[1 -1;-1 1]);    % Kronecker product
Bm=kron(ur.*h_e./6,[2 1;1 2]);

K_midm=Am+k0^2*Bm;

K_diag0m=[K_midm(1:4:end) 0]+[0 K_midm(4:4:end)];
K_diag1m=K_midm(3:4:end);
K_diagM1m=K_midm(2:4:end);
Km=diag(K_diag0m)+diag(K_diag1m,1)+diag(K_diagM1m,-1);

% 1.24
%% -------------------G P-------------------------------------------

%----------------- for E ---------------------------------------
coef_g=a1.*(1./ur-er)*k0^2;
% 
z_gq=Node.'*h_e/2+repmat((z_e1+z_e2)/2,length(Node),1);
g_temp1=(repmat(z_e2,length(Node),1)-z_gq).*exp(-j*k0*z_gq); 
g_temp2=(z_gq-repmat(z_e1,length(Node),1)).*exp(-j*k0*z_gq);
g_e1=Weight*g_temp1.*coef_g/2;  
g_e2=Weight*g_temp2.*coef_g/2;

G=[g_e1 0]+[0 g_e2];
 
P=zeros(size(z_seg));
P(2:end-1)=(1./ur(1:end-1)-1./ur(2:end))*-j*k0*a1.*exp(-j*k0*z_seg(2:end-1));   % inner boundary condition for scattered field

%----------------- for H ---------------------------------------
coef_gm=a1/sqrt(u0/e0).*(1./er-ur)*k0^2;
% 
z_gq=Node.'*h_e/2+repmat((z_e1+z_e2)/2,length(Node),1);
g_temp1m=(repmat(z_e2,length(Node),1)-z_gq).*exp(-j*k0*z_gq); 
g_temp2m=(z_gq-repmat(z_e1,length(Node),1)).*exp(-j*k0*z_gq);
g_e1m=Weight*g_temp1m.*coef_gm/2;  
g_e2m=Weight*g_temp2m.*coef_gm/2;

Gm=[g_e1m 0]+[0 g_e2m];
 
Pm=zeros(size(z_seg));
Pm(2:end-1)=(1./er(1:end-1)-1./er(2:end))*-j*k0*a1/sqrt(u0/e0).*exp(-j*k0*z_seg(2:end-1));   % inner boundary condition for scattered field




%% -----ABC, source at z->-infinity--------------------------
% ---- E---------------
K(1,1)=K(1,1)-j*k0/ur(1);                          % ABC
K(end,end)=K(end,end)-j*k0/ur(end);
% ---- H---------------
Km(1,1)=Km(1,1)-j*k0/er(1);                          % ABC
Km(end,end)=Km(end,end)-j*k0/er(end);

%---------------- result -------------------------------
% ---- E---------------
F_scat=inv(K)*(G+P).';
F_inc=a1*exp(-j*k0*z_seg);
F=F_scat.' +F_inc;
% ---- H---------------
Fm_scat=inv(Km)*(Gm+Pm).';
Fm_inc=a1/sqrt(u0/e0)*exp(-j*k0*z_seg);
Fm=Fm_scat.' +Fm_inc;


%% -------------------------post process----------------------

%----------S paramter------------------

z_in=5e-3;
z_out=45e-3;
[USLS,N_z_in]=min(abs(z_seg-z_in));
[USLS,N_z_out]=min(abs(z_seg-z_out));
S_11=F_scat(N_z_in)/F_inc(N_z_in);
S_21=F(N_z_out)/F_inc(N_z_in);
S_11m=Fm_scat(N_z_in)/Fm_inc(N_z_in);
S_21m=Fm(N_z_out)/Fm_inc(N_z_in);

%---------- Power dissipation----------
%%---------- P_loss calculated by Poynting theorem ---------
    P_e_losse=zeros(size(er));
    P_e_lossm=zeros(size(er));
    
for i=1:length(z_seg)-1
    B_temp=B(:,2*i-1:2*i);
    P_e_losse(i)=w*e0/2*(F(i:i+1)*-imag(B_temp)*(F(i:i+1)'));
    Bm_temp=Bm(:,2*i-1:2*i);
    P_e_lossm(i)=w*u0/2*(Fm(i:i+1)*-imag(Bm_temp)*(Fm(i:i+1)'));
end;
    p_losse=sum(P_e_losse);
    p_lossm=sum(P_e_lossm);
    p_loss=p_losse+p_lossm;
    p_loss=abs(p_loss);
    

    
%
S11e_T(i2)=S_11;
S21e_T(i2)=S_21;
S11m_T(i2)=S_11m;
S21m_T(i2)=S_21m;
P_loss_T(i2)=p_loss;
end
%}
%{
figure(1)
plot(N_swp,(abs(S11e_T)));
hold on;
plot(N_swp,(abs(S21e_T)),'o');
hold on;
plot(N_swp,(abs(S11m_T)),'--');
hold on;
plot(N_swp,(abs(S21m_T)),'.');
legend('S11_e','S21_e','S11_m','S21_m');
xlabel('Segment number N')
ylabel('Absolute value of S parameters')
set(gcf,'Color',[1,1,1]); % White background

figure(2)
plot(N_swp,P_loss_T);
xlabel('Segment number N')
ylabel('total loss')
set(gcf,'Color',[1,1,1]); % White background
%} 
%
figure(1)
plot(f_swp,(abs(S11e_T)),'linewidth', 1.5);
hold on;
plot(f_swp,(abs(S21e_T)));
hold on;
plot(f_swp,abs(S21e_T).^2+abs(S11e_T).^2,'-.');
legend('S11','S21','square sum of S11 and S21');
xlabel('Frequent range/ Hz')
ylabel('Absolute value of S parameters')
set(gcf,'Color',[1,1,1]); % White background
grid on;
figure(2)
plot(f_swp,P_loss_T);
xlabel('Frequent range/ Hz');
ylabel('total loss')
set(gcf,'Color',[1,1,1]); % White background
grid on;
 %}   