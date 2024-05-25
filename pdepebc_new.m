function [pl,ql,pr,qr] = pdepebc_new(xl,ul,xr,ur,t)
    
    dEdx = 1.5;
    zeta = -0.0027;
    mu_water = 1;
    n=0.64;
    e_w = 5;
    c1 = 200;
    D1 = 9.3100*10^-9;
    k_eo2 = (e_w*zeta)*n/(mu_water);
    Jw = -D1*c1 - (c1*(0 + k_eo2)*dEdx);
    pl = ul;
    ql = 0;
    pr = 0;
    qr = 1;
end