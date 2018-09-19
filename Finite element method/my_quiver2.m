function [ output_args ] = my_quiver2( x, y, vx, vy, h)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    N=prod(size(x));
    if nargin < 5
        if(N==1)
            h=1;
        else
            dif_x=diff(reshape(x,[1,prod(size(x))]));
            min_x=min(abs(dif_x(find(dif_x~=0))));
            dif_y=diff(reshape(y,[1,prod(size(y))]));
            min_y=min(abs(dif_y(find(dif_y~=0))));
            h=max([max(abs(vx)) max(abs(vy))])./min([min_x min_y]);
        end
        if nargin < 4
            error(message('my_quiver:NotEnoughInputs'));
      end
    end
    
    
    x=reshape(x,1,N);
    y=reshape(y,1,N);
    vx=reshape(vx,1,N);
    vy=reshape(vy,1,N);
    h=reshape(h,1,prod(size(h)));
       
    R_arrow=1/6;     %arrowhead ratio
    Ax=x+vx./h;
    Ay=y+vy./h;
    Dx=(1-R_arrow)*vx./h+x;
    Dy=(1-R_arrow)*vy./h+y;
    dVx=-0.5774*R_arrow*vy./h;
    dVy=0.5774*R_arrow*vx./h;
    Bx=Dx+dVx;
    By=Dy+dVy;
    Cx=Dx-dVx;
    Cy=Dy-dVy;
    
    line_x=[x;Ax;NaN(size(x))];
    line_y=[y;Ay;NaN(size(y))];
    line_x=reshape(line_x,1,prod(size(line_x)));
    line_y=reshape(line_y,1,prod(size(line_y)));
    VertP_x=[Ax;Bx;Cx];
    VertP_y=[Ay;By;Cy];

    plot(line_x,line_y,'b-');
    patch(VertP_x,VertP_y,'b');
    
end

