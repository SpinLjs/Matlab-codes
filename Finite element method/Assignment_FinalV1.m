clear all;

u0=4*pi*1e-7;   %N/A^2
e0=8.854187817e-12; %F/m 
eta0=sqrt(u0/e0);
beta_my=1/100e-3;

%% -------------- original setting for structure -------
% z=[0 100e-3 200e-3];               % first row: z;  second row: epsilon_r, mu_r
% p_z=[9 1; 1 1];

% z=[0 8e-3 16e-3 34e-3 40e-3];               % first row: z;  second row: epsilon_r, mu_r 
% p_z=[3 1.156 7.89 1.25 ; 1.88 1 7.33 1.15];
% % 
% z=[-3e-3 0 15e-3 22e-3 30e-3 34e-3 37e-3 40e-3 47e-3];
% p_z=[3.97 1 2-0.1j 4.23 6 2.75-0.9j 1 3.97; 3.11 1 2-0.1j 1.11 1 3.2-0.8j 8 3.11 ];

d1=100e-3;
d2=d1/2;
d=d1+d2;
z=[0 d1/2 d1/2+d2 d];               % first row: z;  second row: epsilon_r, mu_r
p_z=[9 1 9; 1 1 1];


%% ----------- pre-processing ----------------------
seg_N=200;
[z_seg, p_z_seg]=ass3_seg(z,p_z,seg_N); % My costom mesh function ass3_seg()
er=p_z_seg(1,1:end);            % meshed epsilon_r
ur=p_z_seg(2,1:end);            % meshed mu_r
h_e=diff(z_seg);                % width of each element
z_e1=z_seg(1:end-1);
z_e2=z_seg(2:end);


%% ---------------------------- A B T D S -------------------------------
%------------------------------- for E ---------------------------
Ae=kron(-1./ur./h_e,[1 -1;-1 1]);    % Kronecker product
Be=kron(er.*h_e./6,[2 1;1 2]);
Te=kron(-1./ur.*h_e./6.*(beta_my.^2),[2 1;1 2]);
De=kron(j*-1./ur.*beta_my,[0 -1;1 0]);
Se=Ae+Te+De;

S_diag0=Se(1:4:end)+[Se(end) Se(4:4:end-4)];
S_diag1=Se(3:4:end-4);
S_diagM1=Se(2:4:end-4);
S=diag(S_diag0)+diag(S_diag1,1)+diag(S_diagM1,-1);
S(1,end)=Se(end-2);
S(end,1)=Se(end-1);

B_diag0=Be(1:4:end)+[Be(end) Be(4:4:end-4)];
B_diag1=Be(3:4:end-4);
B_diagM1=Be(2:4:end-4);
B=diag(B_diag0)+diag(B_diag1,1)+diag(B_diagM1,-1);
B(1,end)=Be(end-2);
B(end,1)=Be(end-1);
%%%%%%%%%%%%%%%%%%%------------here last time -------------------

%------------------------------------for H ------------------------

Ame=kron(-1./er./h_e,[1 -1;-1 1]);    % Kronecker product
Bme=kron(ur.*h_e./6,[2 1;1 2]);
Tme=kron(-1./er.*h_e./6.*(beta_my.^2),[2 1;1 2]);
Dme=kron(j*-1./er.*beta_my,[0 -1;1 0]);
Sme=Ame+Tme+Dme;

Sm_diag0=Sme(1:4:end)+[Sme(end) Sme(4:4:end-4)];
Sm_diag1=Sme(3:4:end-4);
Sm_diagM1=Sme(2:4:end-4);
Sm=diag(Sm_diag0)+diag(Sm_diag1,1)+diag(Sm_diagM1,-1);
Sm(1,end)=Sme(end-2);
Sm(end,1)=Sme(end-1);

Bm_diag0=Bme(1:4:end)+[Bme(end) Bme(4:4:end-4)];
Bm_diag1=Bme(3:4:end-4);
Bm_diagM1=Bme(2:4:end-4);
Bm=diag(Bm_diag0)+diag(Bm_diag1,1)+diag(Bm_diagM1,-1);
Bm(1,end)=Bme(end-2);
Bm(end,1)=Bme(end-1);

%% -------------------G P-------------------------------------------
% all zeros;
%----------------- E H ---------------------------------------
% No excition


%% -----Periodic  Boundary condition--------------------------
% ---- E H---------------
% Periodic boundary conditions, f and its first-order differential equal at
% two ends, order of matrix can be -1. P(zU)+P(zL)=0.

%% ---------------- result -------------------------------
N=5;    % number of eigenvalue (lowest frequencies)
N_shift=0;
% ---- E---------------
[Fe k2_e]=eig(S,B,'vector');
[USLS idx_e]=sort(abs(k2_e));
k2_ed=k2_e(idx_e);
k2_e_s=k2_ed(1+N_shift:N+N_shift);
Fe_p=Fe(:,idx_e(1+N_shift:N+N_shift));
Fe_p=[Fe_p; Fe_p(1,:)];
Fe_final=Fe_p.*exp(-j*beta_my*repmat(z_seg.',1,N));
w=abs(sqrt((-k2_e_s)./u0./e0));

for i=1:N
    figure(i)
    subplot(211);
    plot(z_seg,abs(Fe_final(:,i)));
    plot_vertiline(z);
    xlabel(['z / m ,\omega=' num2str(w(i))]);
    ylabel('Magnitude of E');
    %
    figure(i)
    subplot(212);
    plot(z_seg,phase(Fe_final(:,i)));
    plot_vertiline(z);
    xlabel(['z / m ,\omega=' num2str(w(i))]);
    ylabel('Phase of E');
    %}
end

% ---- H---------------
[Fm k2_m]=eig(Sm,Bm,'vector');
[USLS idx_m]=sort(abs(k2_m));
k2_md=k2_m(idx_m);
k2_m_s=k2_md(1+N_shift:N+N_shift);        % the robin boundary condition turn out to one lowest error eigenvalue(may due to computation error)
Fm_p=Fm(:,idx_m(1+N_shift:N+N_shift));
Fm_p=[Fm_p; Fm_p(1,:)];
Fm_final=Fm_p.*exp(-j*beta_my*repmat(z_seg.',1,N));
w_m=abs(sqrt(-k2_m_s./u0./e0));


for i=1:N
%     figure(i+N)
    figure(i+N)
    subplot(211);
    plot(z_seg,abs(Fm_final(:,i)));
    plot_vertiline(z);
    xlabel(['z / m ,\omega=' num2str(w_m(i))]);
    ylabel('Magnitude of H');
    %
    figure(i+N)
    subplot(212);
    plot(z_seg,phase(Fm_final(:,i)));
    plot_vertiline(z);
    xlabel(['z / m ,\omega=' num2str(w_m(i))]);
    ylabel('Phase of H');
    %}
end


%% ---------- Power dissipation , stored and Q value----------

%% --------analytical solution--------------------------------------------
%{
er1=1;
er2=16;
d1=85e-3;
d2=15e-3;
d=d1+d2;
n1=sqrt(er1);
n2=sqrt(er2);
k0=linspace(0,8/d,400);
beta=1/d*acos(cos(k0.*n1.*d1).*cos(k0.*n2.*d2)-(n1.^2+n2.^2)/(2*n1*n2).*sin(k0*n1*d1).*sin(k0*n2*d2));

figure(3);
plot(k0*d, real(beta)*d);
set(gcf,'Color',[1,1,1]); % White background
figure(1);
plot(real(beta)*d , k0*d);
hold off;
%}
  