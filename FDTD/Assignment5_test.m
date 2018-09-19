% author: Jishen Li
%%--------------------------------------------------------------------------
% physical parameter
c=2.99792458e8;  % m/s
u0=4*pi*1e-7;   %H/m
e0=1/c/c/u0;
%%-------------------------------------------------------------------------
% question setup
lambda_0= 1.3e-6;   %   unit/m;
nb=3.2;     % refraction index
nc=3.4;
n_eff=3.2608;

f=c/lambda_0;   %Hz
betaz=2*pi*f*n_eff*sqrt(e0*u0);

h=0.3e-6;   % /m
WR_l=1.6e-6;
WR_r=1.6e-6;
S= 1e-6;


%%-----------------------------------------------------------------------
% numerical setup
Lx=8e-6; % unit:m
Ly=4e-6;
Lz=100e-6;
dx=2e-8;
dy=2e-8;
dz=5e-7;

%%-----------------------------------------------------------------------
% initialization
x=0:dx:Lx;
y=0:dy:Ly;
z=0:dz:Lz;
E=zeros(length(y),length(x),length(z));  % Ez Hx Hy for update step\

er=nb.^2*ones(length(y),length(x)); % er_distribution
n_xL2=round((Lx/2-S/2)./dx);        % node number for guide boundaries
n_xL1=round((Lx/2-S/2-WR_l)./dx);
n_xR1=round((Lx/2+S/2)./dx);
n_xR2=round((Lx/2+S/2+WR_l)./dx);
n_h1=round((Ly/2-h/2)./dy);
n_h2=round((Ly/2+h/2)./dy);
er(n_h1:n_h2,n_xL1:n_xL2)=nc.^2;    % left waveguide
er(n_h1:n_h2,n_xR1:n_xR2)=nc.^2;  % right waveguide

% C=(2*pi*f).^2*u0*e0.*(er-n_eff^2);
C=(2*pi*f).^2*u0*e0.*er-betaz.^2;   
Ay=1./dz-1i/betaz/dy/dy+1i.*C./4/betaz; % diagonal terms
Bx=1./dz+1i/betaz/dx/dx-1i.*C./4/betaz;
Ax=1./dz-1i/betaz/dx/dx+1i.*C./4/betaz;
By=1./dz+1i/betaz/dy/dy-1i.*C./4/betaz;
Sid_x=1i/2/betaz/dx/dx;     % off-diagonal terms
Sid_y=1i/2/betaz/dy/dy;

%%-------------------------------------------------------------------------------------------------------
% Source initialization

E_max = 1;  % V/m
m_x = Lx/2-S/2-WR_l/2;  % center of gauusian 
m_y = Ly/2;
p=0.405;     % p,q adjust the x,y direction variance
q=0.866;
theta_x=WR_l/2*sqrt(-1/2/log(p));
theta_y=h/2*sqrt(-1/2/log(q));
E_src=E_max.*exp(-1/2*(y.'-m_y).^2/theta_y/theta_y)*exp(-1/2*(x-m_x).^2/theta_x/theta_x);
E(:,:,1)=E_src; % gaussian as source

%%-------------------------------------------------------------------------------------------------------
% update iteration
figure(2);
for n=1:length(z)-1
    if (mod(n,2)==1)    % n+1 step: x--explicit, y--implicit
        for i=2:length(x)-1
            RH1=sum([-Sid_x*ones(length(y)-2,1) Bx(2:end-1,i) -Sid_x*ones(length(y)-2,1)].*E(2:end-1,i-1:i+1,n),2);
%             Tridia1=diag(Ay(2:end-1,i))+diag(Sid_y*ones(1,length(y)-3),-1)+diag(Sid_y*ones(1,length(y)-3),1);
%             E(2:end-1,i,n+1)=inv(Tridia1)*RH1;
            E(2:end-1,i,n+1) = tridiag(Ay(2:end-1,i),Sid_y*ones(1,length(y)),Sid_y*ones(1,length(y)),RH1);
        end
    else
        for j=2:length(y)-1 % n+2 step: x--implicit, y--explicit
            RH2=sum([-Sid_y*ones(1,length(x)-2); By(j,2:end-1); -Sid_y*ones(1,length(x)-2)].*E(j-1:j+1,2:end-1,n),1);
            E(j,2:end-1,n+1) = tridiag(Ax(j,2:end-1),Sid_x*ones(1,length(x)),Sid_x*ones(1,length(x)),RH2);
%             Tridia2=diag(Ax(j,2:end-1))+diag(Sid_x*ones(1,length(x)-3),-1)+diag(Sid_x*ones(1,length(x)-3),1);
%             E(j,2:end-1,n+1)=inv(Tridia2)*RH2.';
        end        
    end
%     E_show=reshape(E((length(y)+1)/2,:,:),length(x),length(z));
%     mesh(z,x,abs(E_show));
%     pause();
end
    

%--------------------------------------------------------------------
% plot results
E_show=reshape(E((length(y)+1)/2,:,:),length(x),length(z));

figure(1)
mesh(x,y,E_src);  
xlabel('x axis / m');
ylabel('y axis / m');
zlabel('E_z magnitude / V/m');
colorbar();
view(2);

figure(2)
mesh(z,x,abs(E_show));
xlabel('x axis / m');
ylabel('z axis / m');
zlabel('E_z magnitude / V/m');
colorbar();
view(2);


