W=2;
d1=W;
d2=2*W;
S=W;
% W=1;
% d1=1;
% d2=9;
% S=4.5;

c=2.99792458e8;  % m/s
u0=4*pi*1e-7;   %H/m
e0=1/c/c/u0;
er=10;

V0=1;      % Voltage on metal
Vg=0;       % ground bonding Voltage
weight=S+W+S;
height=d1+d2;

h_dis=[0.1];    % step h for densest meshes
num_dis=1;

% tao = 1.841:0.0002:1.843;  % over-relaxation factor
tao = 1.84;
n_record = zeros(size(tao));   % iteration number for given structure
n_af_record = zeros(size(tao)); % iteration number for air-filled structure
er_eff=zeros(size(tao));  % effective permittivity
Z0=zeros(size(tao));    % impedance
%%-------------------------------------------------------------------------
% Discretization
h=h_dis(num_dis);    % Diecretization step (unit:1)
x=0:h:weight/2;    % use the 1/2 symmetry
x_total=0:h:weight;
y=0:h:height;
y_total=0:h:height;
LocMetal_y=d1/h+1;          % save the location of metal strip
LocMetal_x=S/h+1:(S+W/2)/h+1;

for n_t = 1:length(tao);
    
    phi=zeros(length(y),length(x));     % the potential phi
    Epsilon=ones(length(y),length(x));  % a matrix Eplison save the relative epsilon for every node
    Q=0;
    
    
    %%------------------------------------------------------------------------
    % Boundary condition & Initialization
    
    phi(LocMetal_y,LocMetal_x)=V0;
    Epsilon(LocMetal_y,LocMetal_x)=NaN;     % use the NaN represnet the metal region
    phi(1:end,1)=Vg;
    Epsilon(1:end,1)=NaN;
    phi(1,1:end)=Vg;
    Epsilon(1,1:end)=NaN;
    phi(end,1:end)=Vg;
    Epsilon(end,1:end)=NaN;
    phi_af=phi;         % Potential of air-filled structure
    Eplison_airfill=Epsilon;    % The air-filled structure need to set er=1 everywhere except on metal
    Epsilon(LocMetal_y,1:end)=0*Epsilon(LocMetal_y,1:end);  % interface
    Epsilon(1:d1/h, 1:end)=er*Epsilon(1:d1/h, 1:end);   % dielectric with er
    
    R=0; % Residue
    
    Q_last=1;
    Q_last_af=1;
    flag1=0;
    flag2=0;
    for n=1:10000
        for i=1:length(y);
            for j=1:length(x);
                if(~isnan(Epsilon(i,j)))    % check in metal or not
                    if(Epsilon(i,j)==0)     % check on interface or not
                        e1t = Epsilon(i+1,j);   % Temporal variable as epsilon 1
                        e2t = Epsilon(i-1,j);   % Temporal variable as epsilon 2
                        R =phi(i,j-1)*(e1t+e2t)+phi(i,j+1)*(e1t+e2t)...
                            +phi(i+1,j)*2*e1t+phi(i-1,j)*2*e2t;
                        R=R/4/(e1t+e2t)-phi(i,j);       % residue expression for Laplacian on interface
                    else if(j~=length(x))
                            R = (phi(i-1,j)+phi(i+1,j)+phi(i,j-1)+phi(i,j+1))/4-phi(i,j);
                            % residue expression for Laplacian on homogeneous node
                        else
                            R = (4*phi(i,j-1)-phi(i,j-2))/3-phi(i,j);
                            % residue expression for Laplacian on symmetry boundary
                        end
                    end
                    phi(i,j)=phi(i,j)+tao(n_t)*R; % SOR expression
                    
                    if(j==length(x))    % same for airfilled structure
                        R = (4*phi_af(i,j-1)-phi_af(i,j-2))/3-phi_af(i,j);
                    else
                        R = (phi_af(i-1,j)+phi_af(i+1,j)+phi_af(i,j-1)+phi_af(i,j+1))/4-phi_af(i,j);
                    end
                    phi_af(i,j)=phi_af(i,j)+tao(n_t)*R;
                end
            end
        end
        Q1=(phi(LocMetal_y,LocMetal_x(1))-phi(LocMetal_y,LocMetal_x(1)-1))*(e1t+e2t)/2;
        Q2=(sum(phi(LocMetal_y,LocMetal_x))-sum(phi(LocMetal_y+1,LocMetal_x)))*e1t;
        Q4=(sum(phi(LocMetal_y,LocMetal_x))-sum(phi(LocMetal_y-1,LocMetal_x)))*e2t;
        Q_mid=(phi(LocMetal_y,LocMetal_x(end))-phi(LocMetal_y+1,LocMetal_x(end)))*e1t...
            +(phi(LocMetal_y,LocMetal_x(end))-phi(LocMetal_y-1,LocMetal_x(end)))*e2t;
        Q=((Q1+Q2+Q4)*2-Q_mid)*e0;  % calculated Q on metal for each iteration
        
        Q1_af=phi_af(LocMetal_y,LocMetal_x(1))-phi_af(LocMetal_y,LocMetal_x(1)-1);
        Q2_af=(sum(phi_af(LocMetal_y,LocMetal_x))-sum(phi_af(LocMetal_y+1,LocMetal_x)));
        Q4_af=(sum(phi_af(LocMetal_y,LocMetal_x))-sum(phi_af(LocMetal_y-1,LocMetal_x)));
        Q_mid_af=phi_af(LocMetal_y,LocMetal_x(end))-phi_af(LocMetal_y+1,LocMetal_x(end))...
            +phi_af(LocMetal_y,LocMetal_x(end))-phi_af(LocMetal_y-1,LocMetal_x(end));
        Q_af=((Q1_af+Q2_af+Q4_af)*2-Q_mid_af)*e0;   % calculated Q on metal of corresponding air-filled structure
        
        Err_q= 1-Q/Q_last;
        Err_q_af= 1-Q_af/Q_last_af;
        
        if(abs(Err_q)<1e-8 & flag1==0)
            n_sl=n;
            flag1=1;
        end
        if(abs(Err_q_af)<1e-8 & flag2==0)
            n_af=n;
            flag2=1;
        end
        if(flag1==1 & flag2==1)
            break;
        end
        
        Q_last=Q;
        Q_last_af=Q_af;
    end
    
    
    phi_12=phi(1:end,end-1:-1:1);       % patch the quarter result to whole area using symmetry
    phi_total=[phi,phi_12];
    phi_af_12=phi_af(1:end,end-1:-1:1);       % patch the quarter result to whole area using symmetry
    phi_af_total=[phi_af,phi_af_12];
    
    figure(1+(num_dis-1)*4);  % 3D with mesh()
    mesh(x_total,y_total,phi_total);
    colorbar();
    figure(2+(num_dis-1)*4);  % contour plot
    contour(x_total,y_total,phi_total);
    colorbar();
    C=Q/V0;     % capacitance for given structure
    
    figure(3+(num_dis-1)*4);  % 3D with mesh()
    mesh(x_total,y_total,phi_af_total);
    colorbar();
    figure(4+(num_dis-1)*4);  % contour plot
    contour(x_total,y_total,phi_af_total);
    colorbar();
    Ca=Q_af/V0;  % capacitance for air-filled structure
    
    n_record(n_t)=n_sl  % iteration number
    n_af_record(n_t)=n_af  % iteration number
    er_eff(n_t)=C/Ca;     % effective permittivity
    Z0(n_t)=sqrt(u0*e0/Ca/C);     % impedance
    
end


