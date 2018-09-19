V0=10;
a=1.5;
b=1;
k_c2d=2; % k_c2d is the course to dense mesh ratio
n_stop=[1000,2000,6000,10000];  % iteration number for converge for each mesh
%%-------------------------------------------------------------------------
% Discretization
h=0.05;    % Diecretization step (unit:m)
x_1=0:h:a;    % distance of adjacent point of discretization is 0.01m
x=x_1(1:(end+1)/2); % Apply quater symmetry
y_1=0:h:b;
y=y_1(1:(end+1)/2);
phi=zeros(length(y),length(x));


%%------------------------------------------------------------------------
% Boundary condition & Initialization
V1=V0;
V2=0;
phi(1,1:end)=V1;
phi(1:end,1)=V2;
phi(2:end,2:end)=(V1+V2)/2;

for num_dis=1:4
    if(num_dis~=1) % from course mesh to dense mesh
        h=h/k_c2d;     % with node distance shorten to 1/k each time;
        x_1=0:h:a;     % update x,y value to dense mesh
        y_1=0:h:b;
        x=x_1(1:(end+1)/2); 
        y=y_1(1:(end+1)/2);
        phi_temp=kron(phi,ones(k_c2d,k_c2d));   % Kron product to change coarse result to dense mesh initial
        phi=phi_temp(k_c2d:end,k_c2d:end);  % discard some extra points
    end
    
for n=1:n_stop(num_dis)
    for i=2:length(y);
        if (i==length(y))
            for j=2:length(x);
                phi(i,j)=(4*phi(i-1,j)-phi(i-2,j))/3;   % Ex=0 at vertical symmetry axis
            end
        else
            for j=2:length(x);
                if (j==length(x))
                    phi(i,j)=(4*phi(i,j-1)-phi(i,j-2))/3;   % Ey=0 at vertical symmetry axis   
                else
                    phi(i,j)=phi(i-1,j)+phi(i+1,j)+phi(i,j-1)+phi(i,j+1); % Laplace equation in difference form
                    phi(i,j)=phi(i,j)/4;
                end
            end
       end
    end
end


end

phi_22=phi(end-1:-1:1,end-1:-1:1);  % patch the quarter result to whole area using symmetry
phi_12=phi(1:end,end-1:-1:1);
phi_21=phi(end-1:-1:1,1:end);
phi_total=[phi,phi_12; phi_21,phi_22];

figure(3);  % 3D with mesh()
mesh(x_1,y_1,phi_total);
colorbar();
figure(4);  % contour plot
contour(x_1,y_1,phi_total);
colorbar();
phi(end,end) % the final value of node in the middle, forv convergence check.
