function [A,B,C] = Draw_my_arrow(ox,oy,vx,vy,h)
%UNTITLED Summary of this function goes here
%   origin(x,y) and vector(x,y) and scale degree h
    
    % Check inputs
    if nargin < 5
      h=1;
      if nargin < 4
        error(message('DrawMyArrow:NotEnoughInputs'));
      end
    end

   R_arrow=1/6;     %arrowhead ratio
   A=[ox+vx./h;oy+vy./h];
   D=[(1-R_arrow)*vx./h+ox;(1-R_arrow)*vy./h+oy];
   dV=[-0.5774*R_arrow*vy./h;0.5774*R_arrow*vx./h];
   B=D+dV;
   C=D-dV;
   hold on;
   plot([ox,A(1)],[oy,A(2)],'b-');
   patch([A(1),B(1) C(1)],[A(2),B(2) C(2)],'b')
end

