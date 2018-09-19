function zd
dx=0.1;
x=1:dx:10;
N=length(x);
y=ones(1,round(N/4));
y(round(N/4)+1:N)=-0.5*ones(1,round(3*N/4));
y1=integral(y,dx);
y2=integral(y1,dx);
subplot 311
h1=plot(x(1),y(1));
axis([1,10,-1,1.5])
subplot 312
h2=plot(x(1),y1(1),'k.');
axis([1,10,-1,3])
subplot 313
h3=plot(x(1),y2(1),'b.');
axis([1,10,0,10])
for ii=2:length(x)
    set(h1,'XData',x(1:ii),'YData',y(1:ii))
    set(h2,'XData',x(1:ii),'YData',y1(1:ii))
    set(h3,'XData',x(1:ii),'YData',y2(1:ii))
    drawnow
    pause(0.01)
end
function [Iy]=integral(y,dx)
total=0;
Iy = y*0;
for i=1:length(y)
    total = total + y(i)*dx;
    Iy(i) = total;
end