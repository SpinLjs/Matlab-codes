%%--------------------------------------------------------------------------
% physical parameter
c=2.99792458e8;  % m/s
u0=4*pi*1e-7;   %H/m
e0=1/c/c/u0;
%%-------------------------------------------------------------------------
% question setup
f=2.5e9;  % HZ
er=4;
ur=1;
R=6e-2; % unit/m
% dis_src=0.6e-2;
dis_src=6e-2;
Ez_strgth=1e5; % V/m
Hx_strgth=Ez_strgth*sqrt(e0/u0);    %A/m
%%-----------------------------------------------------------------------
% numerical setup
Wx=100e-2; % unit:m
Wy=100e-2;
Cntr_Cyl=[Wx/2 Wy/2];
dx=2e-3;
dy=dx;
dt_max=sqrt(e0*u0)/sqrt(1/dx.^2+1/dy.^2);  % unit/s
exp_odr_dt = floor(log10(dt_max));
dt=floor(dt_max/10^exp_odr_dt)*10^exp_odr_dt;

num_PML=10;      % PML thickness
n_e=num_PML;    
n_m=num_PML-1;

val_sigma_e=1e1;  % S/m
val_sigma_m=u0/e0*val_sigma_e; % Omega/m

lambda=c/f;
if (lambda<dx*5)
    disp('spatial discretization step may be too large')
end

%%-----------------------------------------------------------------------
% initialization
x=0:dx:Wx;
y=0:dy:Wy;
Ez=zeros(length(y),length(x));  % Ez Hx Hy for update step\
Hx=Ez;
Hy=Ez;
Ezx=Ez; % Ezx Ezy can be represent with less momery
Ezy=Ez;

er_mtrx=ones(size(Ez)); % er_distribution
for i=1:length(y)
   for j=1:length(x)
       dist_y=Cntr_Cyl(2)-(i-1)*dx;
       dist_x=Cntr_Cyl(1)-(j-1)*dx;
       if(dist_y.^2+dist_x.^2<=R.^2)
           er_mtrx(i,j)=er;
       end
   end
end

% % PML parameter initial
% 
sigma_ex=zeros(size(Ez));
sigma_ey=zeros(size(Ez));
sigma_mx=zeros(size(Ez));
sigma_my=zeros(size(Ez));

% sigma_ex(:,1:n_e+1)=val_sigma_e;
% sigma_ex(:,end-n_e:end)=val_sigma_e;
% sigma_ey(1:n_e+1,:)=val_sigma_e;
% sigma_ey(end-n_e:end,:)=val_sigma_e;
% 
% sigma_mx(:,1:n_e+1)=val_sigma_m;
% sigma_mx(:,end-n_e-1:end)=val_sigma_m;
% sigma_my(1:n_e+1,:)=val_sigma_m;
% sigma_my(end-n_e-1:end,:)=val_sigma_m;

sigma_ex(:,2:n_e+1)=val_sigma_e*ones(length(y),1)*((n_e-(1:n_e))./n_e).^2;
sigma_ex(:,end-n_e:end-1)=val_sigma_e*ones(length(y),1)*((0:n_e-1)./n_e).^2;
sigma_ey(2:n_e+1,:)=((n_e-(1:n_e))./n_e).^2.'*val_sigma_e*ones(1,length(x));
sigma_ey(end-n_e:end-1,:)=((0:n_e-1)./n_e).^2.'*val_sigma_e*ones(1,length(x));

sigma_mx(:,2:n_e+1)=val_sigma_m*ones(length(y),1)*((n_e-(1:n_e))./n_e).^2;
sigma_mx(:,end-n_e:end-1)=val_sigma_m*ones(length(y),1)*((0:n_e-1)./n_e).^2;
sigma_my(2:n_e+1,:)=((n_e-(1:n_e))./n_e).^2.'*val_sigma_m*ones(1,length(x));
sigma_my(end-n_e:end-1,:)=((0:n_e-1)./n_e).^2.'*val_sigma_m*ones(1,length(x));

% % sigma_mx(:,2:n_m+1)=val_sigma_m*ones(length(y)-1,1)*((n_m+1-(1:n_m))./n_m).^2;
% % sigma_mx(:,end-n_m:end-1)=val_sigma_m*ones(length(y)-1,1)*((1:n_m)./n_m).^2;
% % sigma_my(2:n_m+1,:)=((n_m+1-(1:n_m))./n_m).^2.'*val_sigma_m*ones(1,length(x)-1);
% % sigma_my(end-n_m:end-1,:)=((1:n_m)./n_m).^2.'*val_sigma_m*ones(1,length(x)-1);

% Source Initial
idx_src_y=round((Cntr_Cyl(2)-R-dis_src)/dy)+1;
% % Ez(idx_src_y,1:end)=Ez_strgth;
% % Hx(idx_src_y,1:end)=Hx_strgth;
% 
idx_src_x=floor((Cntr_Cyl(2)-R)/dx):ceil((Cntr_Cyl(2)+R)/dx);
% idx_src_x=2+num_PML:size(Ez,1)-num_PML-1;
% Ezx(idx_src_y,idx_src_x)=Ez_strgth/2;
% Ezy(idx_src_y,idx_src_x)=Ez_strgth/2;
% Ez(idx_src_y,idx_src_x)=Ez_strgth;
% % Hx(idx_src_y,idx_src_x(2:end-1))=Hx_strgth;

%%-------------------------------------------------------------------------------------------------------
% update iteration


n_count=0;
for t=dt:dt:2000*dt
    %--------------------------------
    % source
    Ezx(idx_src_y,idx_src_x)=sin(2*pi*f*t)*Ez_strgth/2;
    Ezy(idx_src_y,idx_src_x)=sin(2*pi*f*t)*Ez_strgth/2;
    Ez(idx_src_y,idx_src_x)=Ezx(idx_src_y,idx_src_x)+Ezy(idx_src_y,idx_src_x);
    Ezx(idx_src_y-1,idx_src_x)=sin(2*pi*f*(t-dt))*Ez_strgth/2;
    Ezy(idx_src_y-1,idx_src_x)=sin(2*pi*f*(t-dt))*Ez_strgth/2;
    Ez(idx_src_y-1,idx_src_x)=Ezx(idx_src_y,idx_src_x)+Ezy(idx_src_y,idx_src_x);
    Hx(idx_src_y-1,idx_src_x)=sin(2*pi*f*t-(dt/2))*Hx_strgth;
    Hy(idx_src_y-1,idx_src_x)=0;
%     Ezx(101,101)=sin(2*pi*f*t)*10*Ez_strgth/2;
%     Ezy(101,101)=sin(2*pi*f*t)*10*Ez_strgth/2;
%     Ez(101,101)=sin(2*pi*f*t)*10*Ez_strgth;
    %----------------------------------
    % first 0~dt/2 step

    
    Hx(1:end-1,:)=(2*u0-dt.*sigma_my(1:end-1,:))./(2*u0+dt.*sigma_my(1:end-1,:))...
        .*Hx(1:end-1,:)-2*dt/dy./(2*u0+dt.*sigma_my(1:end-1,:))...
        .*(Ez(2:end,:)-Ez(1:end-1,:));
    
    Hy(:,1:end-1)=(2*u0-dt.*sigma_mx(:,1:end-1))./(2*u0+dt.*sigma_mx(:,1:end-1))...
        .*Hy(:,1:end-1)+2*dt/dx./(2*u0+dt.*sigma_mx(:,1:end-1))...
        .*(Ez(:,2:end)-Ez(:,1:end-1));
    
                
    %----------------------------------
    % second dt/2~dt step
    
    Ezx(2:end-1,2:end-1) = (2*e0.*er_mtrx(2:end-1,2:end-1)-dt.*sigma_ex(2:end-1,2:end-1))./(2*e0.*er_mtrx(2:end-1,2:end-1)+dt.*sigma_ex(2:end-1,2:end-1))...
        .*Ezx(2:end-1,2:end-1) + 2*dt/dx./((2*e0.*er_mtrx(2:end-1,2:end-1)+dt.*sigma_ex(2:end-1,2:end-1)))...
        .*(Hy(2:end-1,2:end-1)-Hy(2:end-1,1:end-2));
    Ezy(2:end-1,2:end-1) = (2*e0.*er_mtrx(2:end-1,2:end-1)-dt.*sigma_ey(2:end-1,2:end-1))./(2*e0.*er_mtrx(2:end-1,2:end-1)+dt.*sigma_ey(2:end-1,2:end-1))...
        .*Ezy(2:end-1,2:end-1) - 2*dt/dy./((2*e0.*er_mtrx(2:end-1,2:end-1)+dt.*sigma_ey(2:end-1,2:end-1)))...
        .*(Hx(2:end-1,2:end-1)-Hx(1:end-2,2:end-1));
    Ez=Ezx+Ezy;

    
    %-------------------------------------------------
    %  plot
     
    n_count=n_count+1;
    if  (mod(n_count,10)==0)
        figure(1);
        mesh(x,y,Ez);
        axis([0 1 0 1 -2e5 2e5]);
%         figure(2);
%         mesh(x,y,Hy); 
%         figure(3);
%         mesh(x,y,Hx);
         pause()
        
    end
end    
figure(1);
mesh(x,y,Ez);
figure(2);
mesh(x,y,Ezx);
figure(3);
mesh(x,y,Ezy);
figure(4);
mesh(x,y,Hy);
figure(5);
mesh(x,y,Hx);
