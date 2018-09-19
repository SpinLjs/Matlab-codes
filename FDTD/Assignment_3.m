W=40e-6; % unit:m
d1=100e-6;
S=20e-6;
m=60e-6;
er_xx=43;
er_yy=28;

c=2.99792458e8;  % m/s
u0=4*pi*1e-7;   %H/m
e0=1/c/c/u0;


V0=1;      % Voltage on metal (unit:V)
Vg=0;       % ground bonding Voltage
weight=2*m+2*S+W;

h_dis=[4 2 1 1/2 1/4]*1e-6;    % step h for densest meshes /um
er_eff=zeros(size(h_dis));  % effective permittivity
Z0=zeros(size(h_dis));    % impedance
for num_dis=1:length(h_dis);
    
    %%-------------------------------------------------------------------------
    % Discretization
    h=h_dis(num_dis);    % Diecretization step (unit:1)
    x=0:h:weight/2;    % use the 1/2 symmetry
    x_total=0:h:weight;
    
    LocReg=[round(m/h+1) round((m+S)/h+1) round(x(end)/h+1)];   % the location of different region
    
    
    %%------------------------------------------------------------------------
    % Initialization
    
    phi_mL=zeros(1,LocReg(1));     % the known potential phi along interface
    phi_W=V0*ones(1,LocReg(3)-LocReg(2));
    rho_SL=zeros(1,LocReg(2)-LocReg(1)); % the known charge density rho along interface
    
    N=length(x);
    k=1:N;
    T=zeros(N,N);
    gamma_k=(k-1/2)/(N+1/2)*pi;     % gamma of k
    lambda_k=abs(2*sin(gamma_k/2));  % sqare root of eigen value
    for i=1:N
        T(i,:)=sqrt(2/(N+1/2))*sin(i.*gamma_k); % eigenvectors to T
    end
    %%---------------------------------------------------------------------
    % anisotropic structure
    C_yx=sqrt(er_yy/er_xx); % the constant added for anisotropic material
    
    A1_diag=C_yx*h./lambda_k.*tanh(lambda_k/C_yx/h*d1); % diagonal elements of A1 and A2
    A2_diag=-C_yx*h./lambda_k*1;
    inv_A1=diag(1./A1_diag);   % inverse of diagonal matrix
    inv_A2=diag(1./A2_diag);
    G=er_yy*T*inv_A1*T'-T*inv_A2*T';    % G
    
    G5=G((LocReg(1)+1):LocReg(2),(LocReg(1)+1):LocReg(2));
    G6=G((LocReg(1)+1):LocReg(2),(LocReg(2)+1):LocReg(3));
    
    phi_SL=-inv(G5)*G6*phi_W';  % calculated potential for the unknown interval
    phi=[phi_mL phi_SL' phi_W];
    rho_mL=e0*G(1:LocReg(1),:)*phi';
    rho_W=e0*G(LocReg(2)+1:LocReg(3),:)*phi';
    rho=[rho_mL' rho_SL rho_W'];
    Q=0;
    for i=1:length(rho_W)-1
         Q=Q+h/2*(rho_W(i)+rho_W(i+1)); % trapezoidal rule for integral
    end
    Q=Q*2;
    C=Q/V0;
    %%------------------------------------------------------------------------
    % corresponding air-filled structure;
    A1_diag_af=h./lambda_k.*tanh(lambda_k/h*d1);
    A2_diag_af=-h./lambda_k*1;
    inv_A1_af=diag(1./A1_diag_af);
    inv_A2_af=diag(1./A2_diag_af);
    G_af=T*inv_A1_af*T'-T*inv_A2_af*T';
    
    G5_af=G_af((LocReg(1)+1):LocReg(2),(LocReg(1)+1):LocReg(2));
    G6_af=G_af((LocReg(1)+1):LocReg(2),(LocReg(2)+1):LocReg(3));
    
    phi_SL_af=-inv(G5_af)*G6_af*phi_W';
    phi_af=[phi_mL phi_SL_af' phi_W];
    rho_mL_af=e0*G_af(1:LocReg(1),:)*phi_af';
    rho_W_af=e0*G_af(LocReg(2)+1:LocReg(3),:)*phi_af';
    rho_af=[rho_mL_af' rho_SL rho_W_af'];
    Q_af=0;
    for i=1:length(rho_W_af)-1
        Q_af=Q_af+h/2*(rho_W_af(i)+rho_W_af(i+1));
    end
    Q_af=Q_af*2;
    Ca=Q_af/V0;
    
    er_eff(num_dis)=C/Ca;     % effective permittivity
    Z0(num_dis)=sqrt(u0*e0/Ca/C);     % impedance
    
    
    phi_total=[phi phi(end-1:-1:1)];  % 1/2 symnmetry
    rho_total=[rho rho(end-1:-1:1)];
    figure(2*num_dis-1) % plot of potential
    plot(x_total,phi_total);
    axis([0 2e-4 -0.2 1.2]);
    xlabel('x axis. unit:m');
    ylabel('\Phi_e.   unit:V');
    
    figure(2*num_dis)   % plot of charge density
    plot(x_total,rho_total);
    axis([0 2e-4 -8e-5 10e-5]);
    xlabel('x axis.  unit:  m');
    ylabel('\rho_s.   unit: C/m');
end

figure()    % Z0 vs. h steps
plot(h_dis,Z0,'-',h_dis,Z0,'o')
set(gca,'XDir','reverse'); 
xlabel('Discretization step h');
ylabel('Z_0 / \Omega');

figure()    % relative permittivity vs. h steps
plot(h_dis,er_eff,'-',h_dis,er_eff,'o');
set(gca,'XDir','reverse'); 
xlabel('Discretization step h');
ylabel('\epsilon_r_e_f_f');
