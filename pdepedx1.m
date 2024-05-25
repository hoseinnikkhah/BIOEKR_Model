function [c,f,s] = pdepedx1(x,t,u,DuDx)
    %u_teta = 1;
    %u_c = 1;
    %u_h = 1;
    %D1 = 2; %Componet i=1 Diffuision
    %Tao = 5; %tortuosity of the porous medium
    R1 = 120; %Mollar rate of componet i=1
    c = 1;
    f = DuDx%((D1/Tao^2)*DuDx) - u(u_teta + u_c + u_h);
    s = R1;
end