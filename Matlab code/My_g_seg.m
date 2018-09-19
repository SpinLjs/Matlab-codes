function [ z_seg] = My_g_seg( z ,seg_N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    z_seg=[];
    diff_z=diff(z);
    
    for i=1:length(z)-1
        z_step=diff_z(i)/seg_N;
        z_temp=z(i):z_step:(z(i+1)-z_step);
        z_seg=[z_seg z_temp];     
    end
    z_seg=[z_seg z(end)];     
end

