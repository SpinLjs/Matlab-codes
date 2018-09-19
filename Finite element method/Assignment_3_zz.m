syms zz  z_e2  z_e1  h_e N_e g_zz;
g_zz=0;
N_e=[(z_e2-zz)./h_e;(zz-z_e1)./h_e];
A_sym=diff(N_e,zz)*(diff(N_e,zz)).';
B_sym=N_e*N_e.';
G_sym=N_e*g_zz;
A=zeros(size(A_sym));
B=zeros(size(B_sym));
G=zeros(size(G_sym));
for(i=1:7)
    A=A+Weight(i)*subs(A_sym,zz,(Node(i)*h_e/2+(z_e1+z_e2)/2));    
    B=B+Weight(i)*subs(B_sym,zz,(Node(i)*h_e/2+(z_e1+z_e2)/2));
    %G=G+Weight(i)*subs(G_sym,zz,(Node(i)*h_e/2+(z_e1+z_e2)/2));   
end
A=subs(A,{z_e1,z_e2,h_e},{z_seg(1:end-1),z_seg(2:end),he_n});
B=subs(B,{z_e1,z_e2,h_e},{z_seg(1:end-1),z_seg(2:end),he_n});


%--------------------------------------------------------------
A_11=1./h_e.^2.*(z_e2-z_e1).*p_z_seg;
A_12=-1./h_e.^2.*(z_e2-z_e1).*p_z_seg;
A_21=A_12;
A_22=1./h_e.^2.*(z_e2-z_e1).*p_z_seg;

z_e1_rep=repmat(z_e1,7,1);
z_e2_rep=repmat(z_e2,7,1);
z_rep=Node'*h_e/2+(z_e1_rep+z_e2_rep)/2;
h_e_rep=repmat(h_e,7,1);

B_11=p_z_seg*Weight*((z_e2_rep-z_rep)./h_e_rep).^2;
B_21=p_z_seg*Weight*(z_e2_rep-z_rep)./h_e_rep.*(z_rep-z_e1_rep)./h_e_rep;
B_12=B_21;
B_22=p_z_seg*Weight*((z_rep-z_e1_rep)./h_e_rep).^2;

K_11=A_11-k0^2*B_11;
K_21=A_21-k0^2*B_21;
K_12=K_21;
K_22=A_22-k0^2*B_22;

K_diag0=[K_11 0]+[0 K_22];
K_diag1=K_12;
K_diagM1=K_21;
K=K_diag0+K_diag1+K_diagM1;