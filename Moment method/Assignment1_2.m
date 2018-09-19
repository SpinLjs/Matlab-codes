clc;
itvx = 1e-2;            % x axis point-to-point distance
itvy = 1e-2;            % y axis point-to-point distance
step=1e-3;              % the width of quadrature strip
f=5e9;                  % frequence
w=2*pi*f;               % angular frequence
u0=4*pi*1e-7;           % Vacuum permeability 
e0=8.854e-12;           % Vacuum permittivity F/m
k0=w*sqrt(e0*u0);

a=6e-1/2;               % given parameters /m
xs=5e-2;
xl=-1.2;
xu=1.2;
nod=[0.56522282050801e-2 0.73403717426523e-1 0.284957404462558 0.619482264084778 0.915758083004698];    % 5 points G-MEW 
weight=[0.210469457918546e-1 0.130705540744447 0.289702301671314 0.350220370120399 0.208324841671986];  % weight   
yl=-a;                  % upper and lower boundary of Strip of Current density
yu=a;
yol=-2*a;               % Observe interval of y axis (need if want to draw 3D picture)
you=2*a;

flag=0;

y=yol:itvy:you;       % need if draw 3D pictures
% y=0;                    % only need y=0 data
x=xl:itvx:xu;
fxy=zeros(length(y),length(x));
for i2=1:length(y)      
for i1=1:length(x)
    I=0;
    sl=yl;              % Ex1 from electric current density
    su=sl+step;
    while(su<yu);       % caluculate integration strip-by-strip 
        if x(i1)==xs    % need to be separete to 2 integrations if encounter singularity
            flag=1;
            xi0=(y(i2)-sl)/(su-sl);     % singularity point
            tao1=xi0*(1-nod);           % get IA
            v1=(su-sl)*tao1+sl
            u1=k0*sqrt((x(i1)-xs).^2+(y(i2)-v1).^2);
            f1=xi0*(su-sl)*besselh(0,2,u1);
            tao2=xi0+(1-xi0)*nod;       % get IB
            v2=(su-sl)*tao2+sl;
            u2=k0*sqrt((x(i1)-xs).^2+(y(i2)-v2).^2);
            f2=(1-xi0)*(su-sl)*besselh(0,2,u2);
            I=I+weight*f1'+weight*f2';        
            
        else
            v=(su-sl)*nod+sl;       % directly use Gauss-MRW without singularity
            u=k0*sqrt((x(i1)-xs).^2+(y(i2)-v).^2);
            f=(su-sl)*besselh(0,2,u);
            I=I+weight*f';
        end
        sl=sl+step;         % go to next integration strip 
        su=su+step;
    end
    fxy(i2,i1)=I;           % get f(x,y);
end
end
nomal_fxy=abs(fxy)/abs(fxy((xs-xl)/itvx+1));    %nomalise f(x,y)
% ph_fxy=phase(fxy);              % get phase;
% figure(1);                      % plot graphs
% title('|f(x,0)|/|f(xs,0)|');
% plot(x,nomal_fxy);
% xlabel('x/m');
% ylabel('|f(x,0)|/|f(xs,0)|');
% figure(2);
% title('Phase of fxy');
% plot(x,ph_fxy);
% xlabel('x/m');
% ylabel('Phase of fxy(in ridians)');

mesh(x,y,nomal_fxy);          % get 3D graphic of magnitude
