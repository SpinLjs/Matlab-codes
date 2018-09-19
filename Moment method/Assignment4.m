% Including customed function MyBesselh/ All parameters' unit are standard
% unit
clear all;
clc;

f=15e9;                 % frequence/ Hz
w=2*pi*f;               % angular frequence
u0=4*pi*1e-7;           % Vacuum permeability 
e0=8.854e-12;           % Vacuum permittivity F/m
k=w*sqrt(e0*u0);        % k0
n_TE=1;         % mode number n of TEn


L=6e-1;         % length of PEC plates   
b=12e-3;        % 1/sqrt(u0*e0)/f/2;  % width between two PEC  
N=600;

k_c=n_TE*pi/b;      % cut-off wavenumbers of TEn
f_c=n_TE/(2*b*sqrt(u0*e0));     % cut-off frequencies
beta=sqrt(k^2-k_c^2);       % propagation constant         
lambda_g=2*pi/beta; % guide wavelength
d=lambda_g/4;

Za=0;       % Za Z-axis beginnig point of PEC
Zb=Za+L;    % Zb Z-axis end point of PEC
Ya=0;       % Ya Y-axis beging value
Yb=Ya+b;    % Yb Y-axis end value
LS=[d b/2];    % Line source locatioj
J0=1;       % Value of J0

zitvob=5e-4;    % step of observe points in z-axis
yitvob=2e-4;    % step of observe points in y-axis    
z=L;  % observe points;
y=-2e-3:yitvob:b+2e-3;
z_ob=repmat(z,length(y),1);
y_ob=repmat(y',1,length(z));

% ro=300;         % use to calculate far-zone field with polar coordinates
% theta=-pi:2*pi/2160:pi;
% [ro2,theta2]=meshgrid(ro,theta);
% z_ob=ro2.*cos(theta2);
% y_ob=ro2.*sin(theta2);
% z=zeros(1,length(ro));
% y=zeros(1,length(theta));

zitv = L/N;
yitv = b/N;

% Gauss-MRW Quadrature(5 points)
nod=[0.56522282050801e-2 0.73403717426523e-1 0.284957404462558 0.619482264084778 0.915758083004698];    % 5 points G-MRW 
weight=[0.210469457918546e-1 0.130705540744447 0.289702301671314 0.350220370120399 0.208324841671986];  % weight

PEC_1=[Za:zitv:Zb-zitv];    % Define each PEC plate with its every segments' pole point coordinate
PEC1=[PEC_1;Ya*ones(1,length(PEC_1));PEC_1+zitv;Ya*ones(1,length(PEC_1))];
PEC_2=[Za:zitv:Zb-zitv];
PEC2=[PEC_2;Yb*ones(1,length(PEC_2));PEC_2+zitv;Yb*ones(1,length(PEC_2))];
PEC_3=[0:yitv:b-yitv];
PEC3=[0*ones(1,length(PEC_3));PEC_3;0*ones(1,length(PEC_3));PEC_3+yitv];
% PEC_4=[0:yitv:b-yitv];
% PEC4=[Zb*ones(1,length(PEC_4));PEC_4;Zb*ones(1,length(PEC_4));PEC_4+yitv];
PEC=[PEC1 PEC2 PEC3];       % The matrix represent all PEC
PEC_ct=[(PEC(1,:)+PEC(3,:))/2;(PEC(2,:)+PEC(4,:))/2];  % center of segements of PEC

Vm=w*u0*J0/4*MyBesselh(PEC_ct(1,:),LS(1),LS(1),PEC_ct(2,:),LS(2),LS(2),k); % Vm vector

zm=repmat(PEC_ct(1,:).',1,length(PEC_ct(1,:))); % wirte as matrix for easily get 
zn_1=repmat(PEC(1,:),length(PEC(1,:)),1);       % matrix terms of Zmn with use mang 
zn_2=repmat(PEC(3,:),length(PEC(3,:)),1);       % for loops
ym=repmat(PEC_ct(2,:).',1,length(PEC(2,:)));
yn_1=repmat(PEC(2,:),length(PEC(2,:)),1);
yn_2=repmat(PEC(4,:),length(PEC(4,:)),1);

Zmn=zeros(length(PEC_ct(1,:)));
Zmn_qrad=zeros([size(Zmn) length(nod)]);      % Zmn matrix at different Xi value 
hn=sqrt((PEC(3,:)-PEC(1,:)).^2+(PEC(4,:)-PEC(2,:)).^2); % hn length of a segment line
hn_Z=repmat(hn,length(hn),1);   
for i=1:length(nod) % Gauss-MRW on Zmn
   Zmn_qrad(:,:,i)=MyBesselh(zm,zn_1,zn_2,ym,yn_1,yn_2,k,nod(i));   % User-defined function Mybessel
   Zmn_qrad(:,:,i)=Zmn_qrad(:,:,i).*(1-diag(ones(1,length(Zmn))));  % off-diagonal terms
   xi0=1/2; % the singularity point of Xi
   % diagonal terms
   Zmn_diag_1=xi0*MyBesselh(PEC_ct(1,:),PEC(1,:),PEC(3,:),PEC_ct(2,:),PEC(2,:),PEC(4,:),k,(1-nod(i))*xi0);  
   Zmn_diag_2=(1-xi0)*MyBesselh(PEC_ct(1,:),PEC(1,:),PEC(3,:),PEC_ct(2,:),PEC(2,:),PEC(4,:),k,xi0+(1-xi0)*nod(i));
   Zmn_qrad(:,:,i)=Zmn_qrad(:,:,i)+diag(Zmn_diag_1+Zmn_diag_2); % total Zmn terms in Gauss-MRW 
   Zmn_qrad(:,:,i)=Zmn_qrad(:,:,i)*weight(i);   % Gauss-MRW weightings 
end
Zmn=sum(Zmn_qrad,3).*hn_Z*(-w*u0/4);    % Get Zmn

I=(inv(Zmn)*Vm.');  % Get I vector

E_inc=-w*u0*J0/4*MyBesselh(z_ob,LS(1),LS(1),y_ob,LS(2),LS(2),k);    % incident field at observe points

Az=zeros([length(y) length(z)]);    % Az megnetic vector potential due to equivalent current density
Azn_temp=zeros([length(y) length(z) length(nod)]);
for n=1:length(I)   % Gauss-MRW on Az
    for i=1:length(nod)
       Azn_temp(:,:,i)=MyBesselh(z_ob,PEC(1,n),PEC(3,n),y_ob,PEC(2,n),PEC(4,n),k,nod(i));
       Azn_temp(:,:,i)=Azn_temp(:,:,i)*weight(i);
    end
    Az=Az+I(n)*sum(Azn_temp,3)*hn(n);
end
Az=Az*(u0/4j);  % Az
E_scat=-j*w*Az; % Scattered field at observe points
E_total=E_inc+E_scat;   % Total field

zy_pl(:,:,1)=repmat(z,length(y),1); % change field in PEC to zero
zy_pl(:,:,2)=repmat(y',1,length(z));
zy_idx=ones(length(y),length(z));

for n=1:length(PEC(1,:))
    index_1=find(abs(zy_pl(:,:,1)-PEC(1,n))<zitv/2 & abs(zy_pl(:,:,2)-PEC(2,n))<yitv/2);
    zy_idx(index_1)=0;
end

E_total=E_total.*zy_idx;

figure(1);  % plot electric fields (inc scat and total) 
mesh(z_ob,y_ob,10*log10(abs(E_total)));
% mesh(z_ob,y_ob,abs(E_total));
figure(2);
mesh(z_ob,y_ob,abs(E_scat));
figure(3);
mesh(z_ob,y_ob,abs(E_inc));

figure(1);  % plot electric fields (inc scat and total) 
plot(y_ob,abs(E_total));
xlabel('z axis/ unit m');
ylabel('Magnitude of the total field')
figure(2);
plot(y_ob,abs(E_scat));
figure(3);
plot(y_ob,abs(E_inc));

figure(1);  % plot electric fields (inc scat and total) 
pcolor(z_ob,y_ob,abs(E_total));
xlabel('z axis/ m');
ylabel('y axis/ m');
shading interp
figure(2);
pcolor(z_ob,y_ob,abs(E_scat));
xlabel('z axis/ m');
ylabel('y axis/ m');
shading interp
figure(3);
pcolor(z_ob,y_ob,abs(E_inc));
xlabel('z axis/ m');
ylabel('y axis/ m');
shading interp


% figure(1);    % plot far-zone field in polar coordinates
% E_dB=10*log10(abs(E_total));
% % E_dB(find(E_dB<0))=0;
% polar(theta.',E_dB);
% figure(2);
% F_pattern=10*log10(abs(E_total./(exp(j*k*ro)./sqrt(ro))));
% % F_pattern(find(F_pattern<0))=0;
% polar(theta.',F_pattern);
% 
% figure(2);
% plot(z_ob,abs(E1),'r');
% hold on;
% plot(z_ob,abs(E2),'g');
% hold on;
% plot(z_ob,abs(E3),'b');
% hold on;
% plot(z_ob,abs(E4),'.r');
% hold on;
% plot(z_ob,abs(E5),'.b');
% xlabel('z axis/ unit m');
% ylabel('Magnitude of the total field')
