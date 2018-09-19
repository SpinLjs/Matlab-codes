function [ output_args ] = my_quiver( x, y, vx, vy, h)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    N=prod(size(x));
    dif_x=diff(reshape(x,[1,prod(size(x))]));
    min_x=min(abs(dif_x(find(dif_x~=0))));
    dif_y=diff(reshape(y,[1,prod(size(y))]));
    min_y=min(abs(dif_y(find(dif_y~=0))));
    if nargin < 5
        if(N==1)
            h=1;
        else
            h=max([max(abs(vx)) max(abs(vy))])./min([min_x min_y]);
        end
        if nargin < 4
            error(message('my_quiver:NotEnoughInputs'));
      end
    end
   
    figure;
    hold on;
    if ((prod(size(h)))==1)
        for i=1:N
            Draw_my_arrow(x(i),y(i),vx(i),vy(i),h);
        end
    else
        for i=1:N
            Draw_my_arrow(x(i),y(i),vx(i),vy(i),h(i));
        end
    end
    hold off;
end

