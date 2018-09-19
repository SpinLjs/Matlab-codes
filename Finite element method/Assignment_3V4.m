clear all;

u0=4*pi*1e-7;   %N/A^2
e0=8.854187817e-12; %F/m 
%% ----- seven points Gauss-legender quaduature -------
Weight=[0.129484966168870 0.279705391489277 0.381830050505119 0.417959183673469 0.381830050505119 0.279705391489277 0.129484966168870];
Node=[-0.949107912342759 -0.741531185599394 -0.405845151377397 0 0.405845151377397 0.741531185599394 0.949107912342759];

a1=exp(j*0);    % Magnitude of input E field

f=40e9;         % frequence
w=2*pi*f;
k0=w*sqrt(u0*e0);

%% -------------- original setting for structure -------
z=[0 20e-3 21e-3 30e-3 31e-3 40e-3];
p_z=[1 1.15 4 1.15 1; 1 1 1 1 1 ];
% z=[0 15e-3 22e-3 30e-3 34e-3 37e-3 40e-3 50e-3];
% p_z=[1 2-0.1j 1 6 2.75-0.9j 1 1; 1 2-0.1j 1 1 3.2-0.8j 8 1 ];
% z=[0 20e-3 30e-3 50e-3];        
% p_z=[1 3-0.1j 1; 1 3-0.1j 1 ];  % first row: epsilon_r, second row: mu_r 

%% ----------- pre-processing ----------------------
seg_N=100;
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

%% ---------------- result -------------------------------
% ---- E---------------
F_scat=inv(K)*(G+P).';
F_inc=a1*exp(-j*k0*z_seg);
F=F_scat.' +F_inc;
% ---- H---------------
Fm_scat=inv(Km)*(Gm+Pm).';
Fm_inc=a1/sqrt(u0/e0)*exp(-j*k0*z_seg);
Fm=Fm_scat.' +Fm_inc;

figure(1)
plot(z_seg,abs(F_scat));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of scattered E')
figure(2)
plot(z_seg,abs(F));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of total E')
figure(3)
plot(z_seg,phase(F_scat));
plot_vertiline(z);
xlabel('z / m')
ylabel('Phase of scattered E')

figure(4)
plot(z_seg,abs(Fm_scat));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of scattered H')
figure(5)
plot(z_seg,abs(Fm));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of total H')

figure(6)
plot(z_seg,phase(Fm_scat));
plot_vertiline(z);
xlabel('z / m')
ylabel('Phase of scattered H')

figure(20)
plot(z_seg,abs(G),'b');
hold on;
plot(z_seg,abs(Gm)*100,'r');
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnetude of G and Gm*100')

figure(21)
plot(z_seg,abs(P),'b');
hold on;
plot(z_seg,abs(Pm)*100,'r');
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnetude of P and Pm*100')
%}

%% -------------------------post process----------------------

%% ----------S paramter------------------

z_in=5e-3;
z_out=45e-3;
[USLS,N_z_in]=min(abs(z_seg-z_in));
[USLS,N_z_out]=min(abs(z_seg-z_out));
S_11=F_scat(N_z_in)/F_inc(N_z_in);
S_21=F(N_z_out)/F_inc(N_z_in);
S_11m=Fm_scat(N_z_in)/Fm_inc(N_z_in);
S_21m=Fm(N_z_out)/Fm_inc(N_z_in);

%{
T(i2)=S_11;
S12_T(i2)=S_21;

plot(f_T,mag2db(abs(T)));
%} 

%% ---------- Power dissipation----------
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
    
    %
    %----------- power loss spatial distribution ----------------
    p_e_loss=P_e_losse+P_e_lossm;
    figure(15);
    plot(z_seg(2:end), real(p_e_loss));
    plot_vertiline(z);
    xlabel('z / m')
    ylabel('Phase of scattered H')
    %}
    
    %{
    %-------------- P_loss calculated with S parameters ------------
    p_loss_S=(1-abs(S_11)^2-abs(S_21)^2);
    P_inc=F_inc(N_z_in)*conj(F_inc(N_z_in))/2/sqrt(u0/e0);
    p_loss_S=p_loss_S*P_inc;
       
    p_loss_Sm=(1-abs(S_11m)^2-abs(S_21m)^2);
    P_incm=Fm_inc(N_z_in)*conj(Fm_inc(N_z_in))/2*sqrt(u0/e0);
    p_loss_Sm=p_loss_Sm*P_incm;
    %}
    
    
%% ----------- equivalent volume and Green's function------------
z_gq=Node.'*h_e/2+repmat((z_e1+z_e2)/2,length(Node),1);
z_ob=min(z):(max(z)-min(z))/1000:max(z);
Ee_sc_ob=zeros(size(z_ob));
Me_sc_ob=zeros(size(z_ob));
Em_sc_ob=zeros(size(z_ob));
Mm_sc_ob=zeros(size(z_ob));

E_inc=a1*exp(-j*k0*z_ob);
H_inc=a1/sqrt(u0/e0)*exp(-j*k0*z_ob);

for i=1:length(z_ob) 
    %--- equivalent electric current caused E an H
    %---- Gee and Ni(expansion function) ----------------
    Psi_ee_t1=(repmat(z_e2,length(Node),1)-z_gq).*exp(-j*k0*abs(z_ob(i)-z_gq)); 
    Psi_ee_t2=(z_gq-repmat(z_e1,length(Node),1)).*exp(-j*k0*abs(z_ob(i)-z_gq));
    %---- Gme  ----------------
    Psi_me_t1=(repmat(z_e2,length(Node),1)-z_gq).*(-sign(z_ob(i)-z_gq).*exp(-j*k0*abs(z_ob(i)-z_gq))); 
    Psi_me_t2=(z_gq-repmat(z_e1,length(Node),1)).*(-sign(z_ob(i)-z_gq).*exp(-j*k0*abs(z_ob(i)-z_gq)));
    %---- Gem  ----------------
    Psi_em_t1=(repmat(z_e2,length(Node),1)-z_gq).*(-sign(z_ob(i)-z_gq).*exp(-j*k0*abs(z_ob(i)-z_gq))); 
    Psi_em_t2=(z_gq-repmat(z_e1,length(Node),1)).*(-sign(z_ob(i)-z_gq).*exp(-j*k0*abs(z_ob(i)-z_gq)));
    %---- Gmm  ----------------
    Psi_mm_t1=(repmat(z_e2,length(Node),1)-z_gq).*exp(-j*k0*abs(z_ob(i)-z_gq)); 
    Psi_mm_t2=(z_gq-repmat(z_e1,length(Node),1)).*exp(-j*k0*abs(z_ob(i)-z_gq));
    
    %---- weights to calculate integral --------
    Psi_ee1=j/2/k0.*(Weight*Psi_ee_t1)/2;
    Psi_ee2=j/2/k0.*(Weight*Psi_ee_t2)/2;
    Psi_me1=1/2.*(Weight*Psi_me_t1)/2;
    Psi_me2=1/2.*(Weight*Psi_me_t2)/2;
    
    Psi_em1=1/2.*(Weight*Psi_em_t1)/2;
    Psi_em2=1/2.*(Weight*Psi_em_t2)/2;
    Psi_mm1=j/2/k0.*(Weight*Psi_mm_t1)/2;
    Psi_mm2=j/2/k0.*(Weight*Psi_mm_t2)/2;
    
    % ----------- some coefficients and combining expansion function with E/H node values-- 
    Ee_seg=(er-1).*(F(1:end-1).*Psi_ee1+F(2:end).*Psi_ee2);
    He_seg=(er-1).*(F(1:end-1).*Psi_me1+F(2:end).*Psi_me2);
    Em_seg=(ur-1).*(Fm(1:end-1).*Psi_em1+Fm(2:end).*Psi_em2);
    Hm_seg=(ur-1).*(Fm(1:end-1).*Psi_mm1+Fm(2:end).*Psi_mm2);
    % --------scattered E H field with green function and equivalent
    % source--------------------------------
    E_sc_ob(i)=j*w*e0*sum(Ee_seg)*j*w*u0+j*w*u0*sum(Em_seg);
    H_sc_ob(i)=j*w*e0*sum(He_seg)+j*w*u0*sum(Hm_seg)*j*w*e0; 
    
end

% --------- total E/H field with Green's function method-----
E_ob=E_inc+E_sc_ob;
H_ob=H_inc+H_sc_ob;

figure(7)
plot(z_ob,abs(E_sc_ob));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of scat E (Green)')
figure(8)
plot(z_ob,abs(E_ob));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of total E (Green)')
figure(9)
plot(z_ob,abs(H_sc_ob));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of scat H (Green)')
figure(10)
plot(z_ob,abs(H_ob));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of total H (Green)')

figure(11)
plot(z_ob,phase(E_sc_ob));
plot_vertiline(z);
xlabel('z / m')
ylabel('pahse of scat E (Green)')
figure(12)
plot(z_ob,phase(H_sc_ob));
plot_vertiline(z);
xlabel('z / m')
ylabel('pahse of scat H (Green)')
%{
figure(14)
subplot(211)
plot(z_ob,abs(E_ob));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of total E (Green)')
subplot(212)
plot(z_ob,abs(H_ob));
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of scat H (Green)')

%}
figure(15)
plot(z_seg,abs(F),'linewidth', 1.5);
hold on;
plot(z_seg,abs(F_scat));
hold on;
plot(z_seg,abs(F_inc),'-.');
plot_vertiline(z);
xlabel('z / m')
ylabel('Magnitude of E fields')
legend('E total','E scattered','E incident' );

figure(16)
plot(z_seg,phase(F),'linewidth', 1.5);
hold on;
plot(z_seg,phase(F_scat));
hold on;
plot(z_seg,phase(F_inc),'-.');
plot_vertiline(z);
xlabel('z / m')
ylabel('Phase of E fields')
legend('E total','E scattered','E incident' );
%}
