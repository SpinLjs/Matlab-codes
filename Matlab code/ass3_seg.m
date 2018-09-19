function [ z_seg, p_z_seg] = ass3_seg( z, p_z ,seg_N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    z_seg=[]; p_z_seg=[];
    dz=diff(z);
    dz_min=min(dz);
    
    for i=1:length(z)-1
        z_step=dz(i)./floor(dz(i)/dz_min*seg_N);
        z_temp=z(i):z_step:(z(i+1)-z_step);
        p_z_temp=p_z(1:end,i)*ones(1,length(z_temp));
        z_seg=[z_seg z_temp];
        p_z_seg=[p_z_seg p_z_temp];        
    end
    
    z_seg=[z_seg z(end)];

end

