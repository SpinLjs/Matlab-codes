clear all;
clc;

% er_real=4.17;                % for convergence check, the number of separated-strips
% loss_tan=0.00;           % loss tangent for Complex permittivity
h=0:2e-5:8e-2;
for i4=1:length(h)
zitv = 2e-5;      % minium interval of z axis 
deltaz=zitv;       % set strip interval equal to step of z(convinient to plot)
zin=0;                  % the observe point Zin
         
er=[38 1 38];   % epsilonr values for three slab separately
% er=[1,er_real(i4),1];  
d=[1e-3 h(i4) 1e-3];           % width for three slabs/unit: m
dlu=[1e-2 1e-2];        % distances between slab and observe points
E0=1;                   % E0 value

u0=4*pi*1e-7;           % Vacuum permeability 
e0=8.854e-12;           % Vacuum permittivity F/m

eta0=sqrt(u0/e0);       % eta0

f=1e10;                  % frequency/Hz
% d_lambda=0:0.01:2.1;
% f=1e9:5e7:1e10;
for i3=1:length(f);     % for question of frequence scanning
   
zout=sum(d)+sum(dlu)+zin;   % the observe point Zin                   
z=zin:zitv:zout;        % Z axis
epsilonr=[er(1)*ones(1,int16(d(1)/deltaz))...      % create a set of relative epsilon point-to-point with z in slab
    er(2)*ones(1,int16(d(2)/deltaz)) er(3)*ones(1,int16(d(3)/deltaz))];

w=2*pi*f(i3);           % angular frequence
k0=w*sqrt(e0*u0);       % k0

                 % next 4 rows are convinient process for easily 
                 % change code to more complex situation
zc_1=z(int16(dlu(1)/zitv+1):int16(dlu(1)/zitv)+int16(sum(d)/zitv)); % specific z axis for slab 
idx=find(epsilonr~=1);  % find the index of slab area of permittivity equal to e0, not take to calculation process. 
zc=zc_1(idx)+deltaz/2;  % Zc(m) without area with permittivity epsilon0.
epslr_zc=epsilonr(idx); % epsilon_r set without terms equal to 1.

I_total=zeros(1,length(zc_1)); % pre-difine scale of I,Zmn,V  
Zmn=zeros(length(zc),length(zc));
V=zeros(length(zc),1);
I=zeros(length(zc),1);
for(m=1:length(zc)) % loops to get elements of Zmn and Vm.
    for(n=1:length(zc))
        if(zc(m)==zc(n))
            Zmn(m,n)=1/(j*w*e0*(epslr_zc(m)-1))-j*eta0/k0*(1-exp(-j*k0*deltaz/2));  %m=n
        else
            Zmn(m,n)=eta0/k0*sin(k0*deltaz/2)*exp(-j*k0*abs(zc(m)-zc(n))); %m>n or m<n
        end
    end
    V(m)=E0*exp(-j*k0*zc(m));   % V vector
end
I=inv(Zmn)*V;   % calculate I vector

z_ob=[z(1) z(end)]; % to get Reflect_coef, Transmission_coef and Power loss, 
                    % we only need Ex at Zin and Zout two points
Ex_scat_ob=zeros(1,2);
for i1=1:length(z_ob);
    Itg1=0;               % For every i1, I1 equal to one strip's integration
    for i2=1:length(I)       % caluculate integration strip-by-strip
        f1=-eta0/k0*sin(k0*deltaz/2)*exp(j*k0*-abs(z_ob(i1)-zc(i2)));
        Itg1=Itg1+I(i2)*f1;   % product of I(i2) with E_scatter_xn;
    end
    Ex_scat_ob(i1)=Itg1; % Scatter field at Zin and Zout.
end

Ex_inc_ob=E0*exp(-j*k0.*z_ob);  % Incident field at Zin and Zout.
Ex_ob=Ex_scat_ob+Ex_inc_ob; % Total field at Zin and Zout.

Reflect_coe(i4)=abs(Ex_scat_ob(1)/Ex_inc_ob(1));    % reflect coefficient
Trans_coe(i4)=abs(Ex_ob(end)/Ex_inc_ob(1));     % transmission coefficient

P_in=1/eta0*abs(Ex_inc_ob(1)).^2;       % enter power
P_out=1/eta0*abs(Ex_ob(end)).^2+1/eta0*abs(Ex_scat_ob(1)).^2; % leaving power
P_loss(i4)=P_in-P_out;  % Power loss

end
end

figure(1)       % Plot required figures
grid on;
% plot(d_lambda,Reflect_coe);
plot(h,10*log10(Reflect_coe));
xlabel('d/\lambda');
ylabel('Reflect coefficient |\Gamma|');
figure(2)
% plot(d_lambda,Trans_coe);
plot(h,10*log10(Trans_coe));
xlabel('d/\lambda');
ylabel('Transmission coefficient|T|');
% figure(3)
% plot(f,P_loss);
% xlabel('d/\lambda');
% ylabel('Power loss');
