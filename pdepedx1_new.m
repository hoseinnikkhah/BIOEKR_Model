function [c,f,s] = pdepedx1_new(x,t,u,DuDx)
    c = n;
    f = D0*DuDx+(u_e + k_eo)*dEdx;
    s = R1;
end