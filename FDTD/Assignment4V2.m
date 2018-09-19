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
Wx=40e-2; % unit:m
Wy=40e-2;
Cntr_Cyl=[Wx/2 Wy/2];
dx=2e-3;
dy=dx;
dt_max=sqrt(e0*u0)/sqrt(1/dx.^2+1/dy.^2);  % unit/s
exp_odr_dt = floor(log10(dt_max));
% dt=floor(dt_max/10^exp_odr_dt)*10^exp_odr_dt;
dt=1e-12;

num_PML=5;      % PML thickness
n_e=num_PML;   
n_m=num_PML-1;

val_sigma_e=1e3;  % S/m
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
% for i=1:length(y)
%    for j=1:length(x)
%        dist_y=Cntr_Cyl(2)-(i-1)*dx;
%        dist_x=Cntr_Cyl(1)-(j-1)*dx;
%        if(dist_y.^2+dist_x.^2<=R.^2)
%            er_mtrx(i,j)=er;
%        end
%    end
% end

% % PML parameter initial
% sigma_ex=zeros(size(Ez));
% sigma_ey=sigma_ex;
% sigma_mx=zeros(size(Hy));
% sigma_my=sigma_mx;

sigma_ex=zeros(size(Ez));
sigma_ey=sigma_ex;
sigma_mx=sigma_ex;
sigma_my=sigma_mx;

sigma_ex(:,2:n_e+1)=val_sigma_e;
sigma_ex(:,end-n_e:end-1)=val_sigma_e;
sigma_ey(2:n_e+1,:)=val_sigma_e;
sigma_ey(end-n_e:end-1,:)=val_sigma_e;

sigma_mx(:,1:n_e)=val_sigma_m;
sigma_mx(:,end-n_e+1:end)=val_sigma_m;
sigma_my(1:n_e,:)=val_sigma_m;
sigma_my(end-n_e+1:end,:)=val_sigma_m;

% sigma_ex(:,2:n_e+1)=val_sigma_e*ones(length(y),1)*((n_e-(1:n_e))./n_e).^2;
% sigma_ex(:,end-n_e:end-1)=val_sigma_e*ones(length(y),1)*((0:n_e-1)./n_e).^2;
% sigma_ey(2:n_e+1,:)=((n_e-(1:n_e))./n_e).^2.'*val_sigma_e*ones(1,length(x));
% sigma_ey(end-n_e:end-1,:)=((0:n_e-1)./n_e).^2.'*val_sigma_e*ones(1,length(x));
% 
% sigma_mx(:,2:n_e+1)=val_sigma_m*ones(length(y),1)*((n_e-(1:n_e))./n_e).^2;
% sigma_mx(:,end-n_e:end-1)=val_sigma_m*ones(length(y),1)*((0:n_e-1)./n_e).^2;
% sigma_my(2:n_e+1,:)=((n_e-(1:n_e))./n_e).^2.'*val_sigma_m*ones(1,length(x));
% sigma_my(end-n_e:end-1,:)=((0:n_e-1)./n_e).^2.'*val_sigma_m*ones(1,length(x));

% % sigma_mx(:,2:n_m+1)=val_sigma_m*ones(length(y)-1,1)*((n_m+1-(1:n_m))./n_m).^2;
% % sigma_mx(:,end-n_m:end-1)=val_sigma_m*ones(length(y)-1,1)*((1:n_m)./n_m).^2;
% % sigma_my(2:n_m+1,:)=((n_m+1-(1:n_m))./n_m).^2.'*val_sigma_m*ones(1,length(x)-1);
% % sigma_my(end-n_m:end-1,:)=((1:n_m)./n_m).^2.'*val_sigma_m*ones(1,length(x)-1);

% Source Initial
idx_src_y=round((Cntr_Cyl(2)-R-dis_src)/dy)+1;
Ez(idx_src_y,2:end-1)=Ez_strgth;
Hx(idx_src_y,2:end-1)=Hx_strgth;

% idx_src_x=floor((Cntr_Cyl(2)-R)/dx):ceil((Cntr_Cyl(2)+R)/dx);
% Ez(idx_src_y,idx_src_x)=Ez_strgth;
% Hx(idx_src_y,idx_src_x(2:end-1))=Hx_strgth;
% Ez(101,101)=Ez_strgth;

%%-------------------------------------------------------------------------------------------------------
% update iteration
Reg_1y=1:n_e;
Reg_2y=n_e+1:size(Ez,1)-n_e-1;
Reg_3y=size(Ez,1)-n_e:size(Ez,1)-1;
Reg_1x=1:n_e;
Reg_2x=n_e+1:size(Ez,2)-n_e-1;
Reg_3x=size(Ez,2)-n_e:size(Ez,2)-1;

Reg_e1y=2:n_e+1;
Reg_e2y=n_e+2:size(Ez,1)-n_e-1;
Reg_e3y=size(Ez,1)-n_e:size(Ez,1)-1;
Reg_e1x=2:n_e+1;
Reg_e2x=n_e+2:size(Ez,2)-n_e-1;
Reg_e3x=size(Ez,2)-n_e:size(Ez,2)-1;

n_count=0;
for t=dt:dt:20*dt
  
    %----------------------------------
    % first 0~dt/2 step
    %  non_PML region       
    Hx(Reg_2y,Reg_2x)=Hx(Reg_2y,Reg_2x)-dt/u0/dy.*(Ez(Reg_2y+1,Reg_2x)-Ez(Reg_2y,Reg_2x));       %(j,i)
    Hy(Reg_2y,Reg_2x)=Hy(Reg_2y,Reg_2x)+dt/u0/dx.*(Ez(Reg_2y, Reg_2x+1)-Ez(Reg_2y,Reg_2x));
    
    
    % PML region
    Hx(Reg_1y,2:end-1)=(2*u0-dt.*sigma_my(Reg_1y,2:end-1))./(2*u0+dt.*sigma_my(Reg_1y,2:end-1))...
        .*Hx(Reg_1y,2:end-1)-2*dt/dy./(2*u0+dt.*sigma_my(Reg_1y,2:end-1))...
        .*(Ez(Reg_1y+1,2:end-1)-Ez(Reg_1y,2:end-1));
    Hx(Reg_3y,2:end-1)=(2*u0-dt.*sigma_my(Reg_3y,2:end-1))./(2*u0+dt.*sigma_my(Reg_3y,2:end-1))...
        .*Hx(Reg_3y,2:end-1)-2*dt/dy./(2*u0+dt.*sigma_my(Reg_3y,2:end-1))...
        .*(Ez(Reg_3y+1,2:end-1)-Ez(Reg_3y,2:end-1));
    Hx(Reg_2y,Reg_1x)=(2*u0-dt.*sigma_my(Reg_2y,Reg_1x))./(2*u0+dt.*sigma_my(Reg_2y,Reg_1x))...
        .*Hx(Reg_2y,Reg_1x)-2*dt/dy./(2*u0+dt.*sigma_my(Reg_2y,Reg_1x))...
        .*(Ez(Reg_2y+1,Reg_1x)-Ez(Reg_2y,Reg_1x));
    Hx(Reg_2y,Reg_3x)=(2*u0-dt.*sigma_my(Reg_2y,Reg_3x))./(2*u0+dt.*sigma_my(Reg_2y,Reg_3x))...
        .*Hx(Reg_2y,Reg_3x)-2*dt/dy./(2*u0+dt.*sigma_my(Reg_2y,Reg_3x))...
        .*(Ez(Reg_2y+1,Reg_3x)-Ez(Reg_2y,Reg_3x));
    
    Hy(Reg_1y,2:end-1)=(2*u0-dt.*sigma_mx(Reg_1y,2:end-1))./(2*u0+dt.*sigma_mx(Reg_1y,2:end-1))...
        .*Hy(Reg_1y,2:end-1)+2*dt/dy./(2*u0+dt.*sigma_mx(Reg_1y,2:end-1))...
        .*(Ez(Reg_1y,3:end)-Ez(Reg_1y,2:end-1));
    Hy(Reg_3y,2:end-1)=(2*u0-dt.*sigma_mx(Reg_3y,2:end-1))./(2*u0+dt.*sigma_mx(Reg_3y,2:end-1))...
        .*Hy(Reg_3y,2:end-1)+2*dt/dy./(2*u0+dt.*sigma_mx(Reg_3y,2:end-1))...
        .*(Ez(Reg_3y,3:end)-Ez(Reg_3y,2:end-1));
    Hy(Reg_2y,Reg_1x)=(2*u0-dt.*sigma_mx(Reg_2y,Reg_1x))./(2*u0+dt.*sigma_mx(Reg_2y,Reg_1x))...
        .*Hy(Reg_2y,Reg_1x)+2*dt/dy./(2*u0+dt.*sigma_mx(Reg_2y,Reg_1x))...
        .*(Ez(Reg_2y,Reg_1x+1)-Ez(Reg_2y,Reg_1x));
    Hy(Reg_2y,Reg_3x)=(2*u0-dt.*sigma_mx(Reg_2y,Reg_3x))./(2*u0+dt.*sigma_mx(Reg_2y,Reg_3x))...
        .*Hy(Reg_2y,Reg_3x)+2*dt/dy./(2*u0+dt.*sigma_mx(Reg_2y,Reg_3x))...
        .*(Ez(Reg_2y,Reg_3x+1)-Ez(Reg_2y,Reg_3x));
    
                
    %----------------------------------
    % second dt/2~dt step
    
    % non_PML region  
    Ez(Reg_e2y,Reg_e2x)=Ez(Reg_e2y,Reg_e2x)+dt/e0/dx./er_mtrx(Reg_e2y,Reg_e2x).*(Hy(Reg_e2y,Reg_e2x)-Hy(Reg_e2y,Reg_e2x-1))...
            -dt/e0/dy./er_mtrx(Reg_e2y,Reg_e2x).*(Hx(Reg_e2y,Reg_e2x)-Hx(Reg_e2y-1,Reg_e2x));
        
    % PML region    
    Ezx(Reg_e1y,2:end-1) = (2*e0-dt.*sigma_ex(Reg_e1y,2:end-1))./(2*e0+dt.*sigma_ex(Reg_e1y,2:end-1))...
        .*Ezx(Reg_e1y,2:end-1) + 2*dt/dx./((2*e0+dt.*sigma_ex(Reg_e1y,2:end-1)))...
        .*(Hy(Reg_e1y,2:end-1)-Hy(Reg_e1y,1:end-2));
    Ezy(Reg_e1y,2:end-1) = (2*e0-dt.*sigma_ey(Reg_e1y,2:end-1))./(2*e0+dt.*sigma_ey(Reg_e1y,2:end-1))...
        .*Ezy(Reg_e1y,2:end-1) - 2*dt/dy./((2*e0+dt.*sigma_ey(Reg_e1y,2:end-1)))...
        .*(Hx(Reg_e1y,2:end-1)-Hx(Reg_e1y-1,2:end-1));
    Ez(Reg_e1y,2:end-1)=Ezx(Reg_e1y,2:end-1)+Ezy(Reg_e1y,2:end-1);
    
    Ezx(Reg_3y,2:end-1) = (2*e0-dt.*sigma_ex(Reg_3y,2:end-1))./(2*e0+dt.*sigma_ex(Reg_3y,2:end-1))...
        .*Ezx(Reg_3y,2:end-1) + 2*dt/dx./((2*e0+dt.*sigma_ex(Reg_3y,2:end-1)))...
        .*(Hy(Reg_3y,2:end-1)-Hy(Reg_3y,1:end-2));
    Ezy(Reg_3y,2:end-1) = (2*e0-dt.*sigma_ey(Reg_3y,2:end-1))./(2*e0+dt.*sigma_ey(Reg_3y,2:end-1))...
        .*Ezy(Reg_3y,2:end-1) - 2*dt/dy./((2*e0+dt.*sigma_ey(Reg_3y,2:end-1)))...
        .*(Hx(Reg_3y,2:end-1)-Hx(Reg_3y-1,2:end-1));
    Ez(Reg_3y,2:end-1)=Ezx(Reg_3y,2:end-1)+Ezy(Reg_3y,2:end-1);
    
    Ezx(Reg_e2y,Reg_e1x) = (2*e0-dt.*sigma_ex(Reg_e2y,Reg_e1x))./(2*e0+dt.*sigma_ex(Reg_e2y,Reg_e1x))...
        .*Ezx(Reg_e2y,Reg_e1x) + 2*dt/dx./(2*e0+dt.*sigma_ex(Reg_e2y,Reg_e1x))...
        .*(Hy(Reg_e2y,Reg_e1x)-Hy(Reg_e2y,Reg_e1x-1));
    Ezy(Reg_e2y,Reg_e1x) = (2*e0-dt.*sigma_ey(Reg_e2y,Reg_e1x))./(2*e0+dt.*sigma_ey(Reg_e2y,Reg_e1x))...
        .*Ezy(Reg_e2y,Reg_e1x) - 2*dt/dy./((2*e0+dt.*sigma_ey(Reg_e2y,Reg_e1x)))...
        .*(Hx(Reg_e2y,Reg_e1x)-Hx(Reg_e2y-1,Reg_e1x));
    Ez(Reg_e2y,Reg_e1x)=Ezx(Reg_e2y,Reg_e1x)+Ezy(Reg_e2y,Reg_e1x);
    
    Ezx(Reg_e2y,Reg_e3x) = (2*e0-dt.*sigma_ex(Reg_e2y,Reg_e3x))./(2*e0+dt.*sigma_ex(Reg_e2y,Reg_e3x))...
        .*Ezx(Reg_e2y,Reg_e3x) + 2*dt/dx./((2*e0+dt.*sigma_ex(Reg_e2y,Reg_e3x)))...
        .*(Hy(Reg_e2y,Reg_e3x)-Hy(Reg_e2y,Reg_e3x-1));
    Ezy(Reg_e2y,Reg_e3x) = (2*e0-dt.*sigma_ey(Reg_e2y,Reg_e3x))./(2*e0+dt.*sigma_ey(Reg_e2y,Reg_e3x))...
        .*Ezy(Reg_e2y,Reg_e3x) - 2*dt/dy./((2*e0+dt.*sigma_ey(Reg_e2y,Reg_e3x)))...
        .*(Hx(Reg_e2y,Reg_e3x)-Hx(Reg_e2y-1,Reg_e3x));
    Ez(Reg_e2y,Reg_e3x)=Ezx(Reg_e2y,Reg_e3x)+Ezy(Reg_e2y,Reg_e3x);
    
    
    %-------------------------------------------------
    %  plot
     
    n_count=n_count+1;
    if  (mod(n_count,1)==0)
        figure(1);
        mesh(x,y,Ez);
        figure(2);
        mesh(x,y,Hy);
        figure(3);
        mesh(x,y,Hx);
        pause()
        
    end
end    
figure(1);
mesh(x,y,Ez);
figure(2);
mesh(x,y,Hy);
figure(3);
mesh(x,y,Hx);
   