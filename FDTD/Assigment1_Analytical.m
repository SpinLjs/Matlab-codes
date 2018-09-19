
% Analytical Solution
%%--------------------------------------------------------------------------
a=1.5;
b=1;
x=0:0.01:a;    % distance of adjacent point of discretization is 0.01m
y=0:0.01:b;
V0=10;
phi=zeros(length(y),length(x));
for n=1:100  % The series order n is limited to 200;
    if(1-(-1)^n==2)
        p1=sinh(n*pi*y./a);
        p2=sinh(n*pi/a*(y-b));
        p3=(p1-p2)'*sin(n*pi*x/a);
        phi=phi+4*V0/n/pi/sinh(n*pi*b./a)*p3;
    end
end
figure(1);
mesh(x,y,phi);
colorbar();
figure(2);
contour(x,y,phi);
colorbar();
phi((length(y)+1)/2,(length(x)+1)/2)