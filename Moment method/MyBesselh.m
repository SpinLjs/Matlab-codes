function [ output_args ] = MyBesselh(zm,zn_1,zn_2,ym,yn_1,yn_2,k,xi)
% User-define function for easily calculate Hankel function in form of (1-14)
%   Detailed explanation goes here
thrd=1e-10; % ignoring float point number error lower than 1e-10;
if nargin==7    % default set Xi=0
    xi=0;
end
A=k*sqrt((zm-zn_1-xi*(zn_2-zn_1)).^2+(ym-yn_1-xi*(yn_2-yn_1)).^2);  % parameter in (1-14)
A=floor(A./thrd)*thrd;
output_args=besselh(0,2,A); % Hankel function(0-order 2nd-class)
end

