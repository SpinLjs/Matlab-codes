clear all;

u0=4*pi*1e-7;   %N/A^2
e0=8.854187817e-12; %F/m 
eta0=sqrt(u0/e0);

%% ----- seven points Gauss-legender quaduature -------
Weight=[0.129484966168870 0.279705391489277 0.381830050505119 0.417959183673469 0.381830050505119 0.279705391489277 0.129484966168870];
Node=[-0.949107912342759 -0.741531185599394 -0.405845151377397 0 0.405845151377397 0.741531185599394 0.949107912342759];


%% -------------- original setting for structure -------
z=[0 4e-3 6e-3];               % first row: z;  second row: epsilon_r, mu_r
p_z=[2.5 9; 1 1];

% z=[0 8e-3 16e-3 34e-3 40e-3];               % first row: z;  second row: epsilon_r, mu_r 
% p_z=[3 1.156 7.89 1.25 ; 1.88 1 7.33 1.15];
% % 
% z=[0 15e-3 22e-3 30e-3 34e-3 37e-3 40e-3 50e-3];
% p_z=[1.48 2-0.1j 1 6 2.75-0.9j 1 1; 1 2-0.1j 1 1 3.2-0.8j 8 1 ];

%% ----------- pre-processing ----------------------
seg_N=600;
[z_seg, p_z_seg]=ass3_seg(z,p_z,seg_N); % My costom mesh function ass3_seg()
er=p_z_seg(1,1:end);            % meshed epsilon_r
ur=p_z_seg(2,1:end);            % meshed mu_r
h_e=diff(z_seg);                % width of each element
z_e1=z_seg(1:end-1);
z_e2=z_seg(2:end);


%% ---------------------------- A B -------------------------------
%------------------------------- for E ---------------------------
Ae=kron(-1./ur./h_e,[1 -1;-1 1]);    % Kronecker product
Be=kron(er.*h_e./6,[2 1;1 2]);

A_diag0=[Ae(1:4:end) 0]+[0 Ae(4:4:end)];
A_diag1=Ae(3:4:end);
A_diagM1=Ae(2:4:end);
A=diag(A_diag0)+diag(A_diag1,1)+diag(A_diagM1,-1);
A=A(2:end-1,2:end-1);

B_diag0=[Be(1:4:end) 0]+[0 Be(4:4:end)];
B_diag1=Be(3:4:end);
B_diagM1=Be(2:4:end);
B=diag(B_diag0)+diag(B_diag1,1)+diag(B_diagM1,-1);
B=B(2:end-1,2:end-1);

%------------------------------------for H ------------------------

Ame=kron(-1./er./h_e,[1 -1;-1 1]);    % Kronecker product
Bme=kron(ur.*h_e./6,[2 1;1 2]);

Am_diag0=[Ame(1:4:end) 0]+[0 Ame(4:4:end)];
Am_diag1=Ame(3:4:end);
Am_diagM1=Ame(2:4:end);
Am=diag(Am_diag0)+diag(Am_diag1,1)+diag(Am_diagM1,-1);

Bm_diag0=[Bme(1:4:end) 0]+[0 Bme(4:4:end)];
Bm_diag1=Bme(3:4:end);
Bm_diagM1=Bme(2:4:end);
Bm=diag(Bm_diag0)+diag(Bm_diag1,1)+diag(Bm_diagM1,-1);

%% -------------------G P-------------------------------------------
% all zeros;
%----------------- for E ---------------------------------------
%----------------- for H ---------------------------------------


%% -----PEC Boundary condition--------------------------
% ---- E---------------
% Dirichlet BC, apply to A and B;
% ---- H---------------
% Robin BC, turn out to be zeros;

%% ---------------- result -------------------------------
N=4;    % number of eigenvalue (lowest frequencies) 
% ---- E---------------
[Fe k2_e]=eig(A,B,'vector');
[USLS idx_e]=sort(abs(k2_e));
k2_ed=k2_e(idx_e);
k2_e_s=k2_ed(1:N);
Fe_show=[zeros(1,N); Fe(:,idx_e(1:N)); zeros(1,N)];
w=abs(sqrt((-k2_e_s)./u0./e0));

for i=1:N
    figure(i)
    subplot(211);
    plot(z_seg,abs(Fe_show(:,i)));
    plot_vertiline(z);
    xlabel(['z / m ,\omega=' num2str(w(i))]);
    ylabel('Magnitude of total E');
    %
    figure(i+N)
    subplot(211);
    plot(z_seg,phase(Fe_show(:,i)));
    plot_vertiline(z);
    xlabel(['z / m ,\omega=' num2str(w(i))]);
    ylabel('Phase of total E');
    %}
end
% ---- H---------------
[Fm k2_m]=eig(Am,Bm,'vector');
[USLS idx_m]=sort(abs(k2_m));
k2_md=k2_m(idx_m);
k2_m_s=k2_md(2:N+1);        % the robin boundary condition turn out to one lowest error eigenvalue(may due to computation error)
Fm_show=Fm(:,idx_m(2:N+1));
w_m=abs(sqrt(-k2_m_s./u0./e0));

for i=1:N
%     figure(i+N)
    figure(i)
    subplot(212);
    plot(z_seg,abs(Fm_show(:,i)));
    plot_vertiline(z);
    xlabel(['z / m ,\omega=' num2str(w_m(i))]);
    ylabel('Magnitude of H');
    %
    figure(i+N)
    subplot(212);
    plot(z_seg,phase(Fm_show(:,i)));
    plot_vertiline(z);
    xlabel(['z / m ,\omega=' num2str(w_m(i))]);
    ylabel('Phase of H');
    %}
end


%% ---------- Power dissipation , stored and Q value----------
%
%%---------- P_loss calculated by Poynting theorem ---------
Q_value=zeros(1,N);
for i2=1:N  
    E=Fe_show(:,i2).';
    Hm=Fm_show(:,i2).';
    W_e_losse=zeros(size(er));
    W_e_lossm=zeros(size(er));
    W_e_storee=zeros(size(er));
    W_e_storem=zeros(size(er));

    rate_em=eta0;
    
    for i=1:length(z_seg)-1
        B_temp=Be(:,2*i-1:2*i);
        W_e_losse(i)=e0/2*(E(i:i+1)*-imag(B_temp)*(E(i:i+1)'));
        W_e_storee(i)=e0/4*(E(i:i+1)*real(B_temp)*(E(i:i+1)'));
        Bm_temp=Bme(:,2*i-1:2*i);
        W_e_lossm(i)=u0/2*((Hm(i:i+1)./rate_em)*-imag(Bm_temp)*(Hm(i:i+1)'./rate_em));
        W_e_storem(i)=u0/4*((Hm(i:i+1)./rate_em)*real(Bm_temp)*(Hm(i:i+1)'./rate_em));
    end;
    W_losse=sum(W_e_losse);
    W_lossm=sum(W_e_lossm);
    W_storee=sum(W_e_storee);
    W_storem=sum(W_e_storem);
    W_loss=W_losse+W_lossm;
    W_store=W_storee+W_storem;
    
    Q_value(i2)=real(W_store./W_loss);    
end
    %}
  Q_value;


%% --------analytical solution--------------------------------------------
%{
er1=2.5;
er2=9;
d1=4e-3;
d2=2e-3;
n_w=1;
tao=sqrt(1*u0/er1/e0)*tan(w(n_w)*sqrt(1*u0*er1*e0)*d1)+sqrt(1*u0/er2/e0)*tan(w(n_w)*sqrt(1*u0*er2*e0)*d2)
k1=w(n_w)*sqrt(e0*er1*u0*1);
k2=w(n_w)*sqrt(e0*er2*u0*1);
z1=0:2e-5:d1-2e-5;
z2=d1:2e-5:d1+d2;
A1=1;
E_anal=[A1*sin(k1.*z1) A1*sin(k1*d1)/sin(k2*d2)*sin(k2.*(d1+d2-z2))];
figure(10);
plot([z1 z2], abs(E_anal));
set(gcf,'Color',[1,1,1]); % White background
hold on;
%}