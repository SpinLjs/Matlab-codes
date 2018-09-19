wV0=10;
a=1.5;
b=1;

%%------------------------------------------------------------------------
% Discretization
h=0.01;    % Diecretization step 0.05m
x=0:h:a;    % distance of adjacent point of discretization is 0.01m
y=0:h:b;
phi=zeros(length(y),length(x));


%%------------------------------------------------------------------------
% Boundary condition & Initialization
V1=V0;
V2=0;
V3=V0;
V4=0;
phi(1,1:end)=V1;
phi(end,1:end)=V3;
phi(1:end,1)=V4;
phi(1:end,end)=V2;
phi(2:end-1,2:end-1)=(V1+V2+V3+V4)/4;

%laplace_op1=[0 1 0;1 -4 1;0 1 0];
%laplace_op2=[1 4 1;4 -20 4 1; 1 4 1];
for n=1:5000
    for i=2:length(y)-1;
       for j=2:length(x)-1;
           phi(i,j)=phi(i-1,j)+phi(i+1,j)+phi(i,j-1)+phi(i,j+1);
           phi(i,j)=phi(i,j)/4;
       end
    end
end

figure(3);
mesh(x,y,phi);
colorbar();
figure(4);
contour(x,y,phi);
colorbar();