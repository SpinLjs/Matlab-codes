clear all;

u0=4*pi*1e-7;   %N/A^2
e0=8.854187817e-12; %F/m 
eta0=sqrt(u0/e0);


%% -------------- original setting for structure -------
% 
d1=100e-3;
d2=d1/2;
d=d1+d2;
z=[0 d1/2 d1/2+d2 d];               % first row: z;  second row: epsilon_r, mu_r
p_z=[9 1 9; 1 1 1];

% % 
% z=[-3e-3 0 15e-3 22e-3 30e-3 34e-3 37e-3 40e-3 47e-3];
% p_z=[3.97 1 2-0.1j 4.23 6 2.75-0.9j 1 3.97; 3.11 1 2-0.1j 1.11 1 3.2-0.8j 8 3.11 ];
% d=50e-3;

%% ----------- pre-processing ----------------------
seg_N=100;

[z_seg, p_z_seg]=ass3_seg(z,p_z,seg_N); % My costom mesh function ass3_seg()
er=p_z_seg(1,1:end);            % meshed epsilon_r
ur=p_z_seg(2,1:end);            % meshed mu_r
h_e=diff(z_seg);                % width of each element
z_e1=z_seg(1:end-1);
z_e2=z_seg(2:end);
%% -----------------------sweep---------------------------

N_beta=100;
beta_set=linspace(0,pi,N_beta)/d;

beta_k_e=zeros(10,N_beta)*NaN;
beta_k_m=zeros(10,N_beta)*NaN;
bk_e_line=[];
bk_m_line=[];

for i2=1:1:length(beta_set);
beta_my=beta_set(i2);
%% ---------------------------- A B T D S-------------------------------
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
k_limit=25/d;
% ---- E---------------


[Fe k2_e]=eig(S,B,'vector');
% idx_esort=find(real(k2_e)<0);
% k2_e_sort=k2_e(idx_esort);
% k2_e_sort=k2_e;
[USLS idx_e]=sort(abs(k2_e));
k2_ed=k2_e(idx_e);
k0_e=sqrt(-k2_ed);
k0_e_s=k0_e(find(k0_e<=k_limit));

% k0_e_s=abs(k0_e_s);
beta_k_e(1:length(k0_e_s),i2)=k0_e_s;
bk_e_line=[bk_e_line [d*(k0_e_s.');beta_my*d*ones(1,length(k0_e_s))]];


% ---- H---------------
[Fm k2_m]=eig(Sm,Bm,'vector');
% idx_msort=find(real(k2_m)<0);
% k2_m_sort=k2_m(idx_msort);
% k2_m_sort=k2_m;
[USLS idx_m]=sort(abs(k2_m));
k2_md=k2_m(idx_m);
k0_m=sqrt(-k2_md);
k0_m_s=k0_m(find(k0_m<=k_limit)); 

% k0_m_s=abs(k0_m_s);
beta_k_m(1:length(k0_m_s),i2)=k0_m_s;
bk_m_line=[bk_m_line [d*(k0_m_s).';beta_my*d*ones(1,length(k0_m_s))]];

end

figure(1)
set(gcf,'Color',[1,1,1]); % White background
plot(bk_e_line(2,:),bk_e_line(1,:),'bx');
xlabel('Value of \beta*d');
ylabel('k_0*d for E');
hold on;


%% --------------------- show band ----------------------------------
%{
figure(3)
set(gcf,'Color',[1,1,1]); % White background
xlabel('Value of \beta*d');
ylabel('k_0*d for E');
hold on;
for i=1:length(bk_e_line(1,:))
    plot_horizonline(bk_e_line(1,i));
end
hold off;
%}
%%
%
figure(2)
set(gcf,'Color',[1,1,1]); % White background
plot(bk_m_line(2,:),bk_m_line(1,:),'bx');
xlabel('Value of \beta*d');
ylabel('k_0*d for H');
hold on;
%}

%% --------------------- show band ----------------------------------
%{
for i=1:length(bk_m_line(1,:))
    plot_horizonline(bk_m_line(1,i));
end
hold off;
%}

%% --------analytical solution--------------------------------------------
%
er1=9;
er2=1;
d1=100e-3;
d2=d1/2;
d=d1+d2;
n1=sqrt(er1);
n2=sqrt(er2);
k0=linspace(0,25/d,5000);
beta=1/d*acos(cos(k0.*n1.*d1).*cos(k0.*n2.*d2)-(n1.^2+n2.^2)/(2*n1*n2).*sin(k0*n1*d1).*sin(k0*n2*d2));

idx_cplx=find(imag(beta)~=0)
beta(idx_cplx)=NaN;
    
figure(3);
plot(k0*d, real(beta)*d);
set(gcf,'Color',[1,1,1]); % White background
xlabel('k_0*d');
ylabel('\beta*d');
figure(1);
plot(real(beta)*d , k0*d);
hold off;
figure(2);
plot(real(beta)*d , k0*d);
hold off;
%}
