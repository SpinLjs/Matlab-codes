W=2;
d1=W;
d2=2*W;
S=W;
c=2.99792458e8;  % m/s
u0=4*pi*1e-7;   %H/m
e0=1/c/c/u0;
er=10;

V0=10;  % Voltage on metal strip
Vg=0;   % Voltage on bounding frame
weight=S+W+S;
height=d1+d2;

% h_dis=[0.1];  % step h for increasing dense meshes
h_dis=[1 0.5 0.25 0.125];    % step h for increasing dense meshes
n_record = zeros(size(h_dis));   % iteration number
er_eff=zeros(size(h_dis));  % effective permittivity 
Z0=zeros(size(h_dis));    % impedance
%%-------------------------------------------------------------------------
% Discretization
for num_dis=1:length(h_dis)
    h=h_dis(num_dis);    % Diecretization step (unit:1)
    x=0:h:weight;    % distance of adjacent point of discretization
    y=0:h:height;
    phi=zeros(length(y),length(x));     % the potential phi
    Epsilon=ones(length(y),length(x));  % a matrix Eplison save the relative epsilon for every node
    Q=0;
    LocMetal_y=d1/h+1;  % save the location of metal strip 
    LocMetal_x=S/h+1:(S+W)/h+1;
    
    %%------------------------------------------------------------------------
    % Boundary condition & Initialization
    
    phi(LocMetal_y,LocMetal_x)=V0;  
    Epsilon(LocMetal_y,LocMetal_x)=NaN;     % use the NaN represnet the metal region
    phi(1:end,1)=Vg;    
    Epsilon(1:end,1)=NaN;
    phi(1:end,end)=Vg;
    Epsilon(1:end,end)=NaN;
    phi(1,1:end)=Vg;
    Epsilon(1,1:end)=NaN;
    phi(end,1:end)=Vg;
    Epsilon(end,1:end)=NaN;
    phi_af=phi;         % Potential of air-filled structure
    Eplison_airfill=Epsilon;    % The air-filled structure need to set er=1 everywhere except on metal
    Epsilon(LocMetal_y,1:end)=0*Epsilon(LocMetal_y,1:end);  % interface
    for i= 2:d1/h
        Epsilon(i, 2:end-1)=er;         % dielectric with er
    end
    
   
    J_temp = phi;   % Matrix saving last interation result in Jacobian Method
    Jt_af = phi_af;
    
    Q_last=1; 
    Q_last_af=1;
    for n=1:100000   % the upper limit of iteration set to be 10000
        for i=1:length(y);
            for j=1:length(x);  % check in metal or not
                if(~isnan(Epsilon(i,j)))    % check on metal or not   
                    if(Epsilon(i,j)~=0)     % check on interface or not
                        phi(i,j)=J_temp(i-1,j)+J_temp(i+1,j)+J_temp(i,j-1)+J_temp(i,j+1); % Laplacian on homogeneous node
                        phi(i,j)=phi(i,j)/4;
                    else
                        e1t = Epsilon(i+1,j);   % Temporal variable as epsilon 1
                        e2t = Epsilon(i-1,j);   % Temporal variable as epsilon 2
                        phi(i,j)=J_temp(i,j-1)*(e1t+e2t)+J_temp(i,j+1)*(e1t+e2t)...
                            +J_temp(i+1,j)*2*e1t+J_temp(i-1,j)*2*e2t;   % Laplacian on interface
                        phi(i,j)=phi(i,j)/4/(e1t+e2t);
                    end
                    phi_af(i,j)=Jt_af(i-1,j)+Jt_af(i+1,j)+Jt_af(i,j-1)+Jt_af(i,j+1); % Same for the air-filled structure
                    phi_af(i,j)=phi_af(i,j)/4;
                end
            end
        end
        
        J_temp=phi;
        Jt_af=phi_af;
        
        Q1=(phi(LocMetal_y,LocMetal_x(1))-phi(LocMetal_y,LocMetal_x(1)-1))*(e1t+e2t)/2;
        Q3=(phi(LocMetal_y,LocMetal_x(end))-phi(LocMetal_y,LocMetal_x(end)+1))*(e1t+e2t)/2;
        Q2=(sum(phi(LocMetal_y,LocMetal_x))-sum(phi(LocMetal_y+1,LocMetal_x)))*e1t;
        Q4=(sum(phi(LocMetal_y,LocMetal_x))-sum(phi(LocMetal_y-1,LocMetal_x)))*e2t;
        Q=(Q1+Q2+Q3+Q4)*e0; % calculated Q on metal for each iteration 
        
        Q1_af=phi_af(LocMetal_y,LocMetal_x(1))-phi_af(LocMetal_y,LocMetal_x(1)-1);
        Q3_af=-(phi_af(LocMetal_y,LocMetal_x(end)+1)-phi_af(LocMetal_y,LocMetal_x(end)));
        Q2_af=(sum(phi_af(LocMetal_y,LocMetal_x))-sum(phi_af(LocMetal_y+1,LocMetal_x)));
        Q4_af=(sum(phi_af(LocMetal_y,LocMetal_x))-sum(phi_af(LocMetal_y-1,LocMetal_x)));
        Q_af=(Q1_af+Q2_af+Q3_af+Q4_af)*e0; % calculated Q on metal of corresponding air-filled structure
        Err_q= 1-Q/Q_last;  % use the Q change ratio in two iterations to indicate convergency
        Err_q_af= 1-Q_af/Q_last_af;
        if(abs(Err_q)<1e-8 & abs(Err_q_af)<1e-8)    % if the ratio <1e-8, regard it as converged
            break;
        end
        Q_last=Q;
        Q_last_af=Q_af;
    end
    
    figure(1+(num_dis-1)*4);  % 3D with mesh()
    mesh(x,y,phi);
    colorbar();
    figure(2+(num_dis-1)*4);  % contour plot
    contour(x,y,phi);
    colorbar();
    C=Q/V0;     % capacitance for given structure
    
    figure(3+(num_dis-1)*4);  % 3D with mesh()
    mesh(x,y,phi_af);
    colorbar();
    figure(4+(num_dis-1)*4);  % contour plot
    contour(x,y,phi_af);
    colorbar();
    Ca=Q_af/V0; % capacitance for air-filled structure
    
    n_record(num_dis)=n  % iteration number
    er_eff(num_dis)=C/Ca     % effective permittivity 
    Z0(num_dis)=sqrt(u0*e0/Ca/C)     % impedance
    %}
end

% Richardson extrapolation result
er_eff_RE =((h_dis(end-1)/h_dis(end)).^2*er_eff(end)-er_eff(end-1))/((h_dis(end-1)/h_dis(end)).^2-1)
Z0_RE = ((h_dis(end-1)/h_dis(end)).^2*Z0(end)-Z0(end-1))/((h_dis(end-1)/h_dis(end)).^2-1)

figure()    % show convergence
plot(h_dis,Z0,'-',h_dis,Z0,'o');
axis([0 1 30 45]);
set(gca,'XDir','reverse'); 
xlabel('Discretization step h');
ylabel('Z_0 / \Omega');

figure()
plot(h_dis,er_eff,'-',h_dis,er_eff,'o');
axis([0 1 5 6]);
set(gca,'XDir','reverse'); 
xlabel('Discretization step h');
ylabel('\epsilon_r_e_f_f');

