clear all;
clc;
z=1:1:600;
Jx=cos(pi/120*(z-260)).*[zeros(1,200) ones(1,120) zeros(1,280)]; %
My=376.7*[zeros(1,400) ones(1,80) zeros(1,120)];
figure(1);
plot(z,Jx);
title('electrical current density')
xlabel('z/mm');
ylabel('electrical current density');
figure(2);
plot(z,My);
title('magnetic current density')
xlabel('z/mm');
ylabel('magnetic current density');
