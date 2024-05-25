clc;
clear;
close all;

%coefficents
F = 9.64853321233100184*10^4; %Faraday Constant [C/mol]
z1 = 5; %charge number of i=1 componet
R1 = 120; %Mollar rate of componet i=1
R = 8.314; %Gas constant [J/mol.k-1]
T = 300; %Temp in Kelvin
D1 = 2; %Componet i=1 Diffuision
Tao = 5; %tortuosity of the porous medium
v1 = D1/(R*T); %Mobility of componet i=1
u_teta = 1;
u_c = 1;
u_h = 1;


%defining pde charachterists
m = 0;
dx = 0.1 ;
dt = 0.1 ;
x = 0 : dx : 1 ;
t = 0 : dt : 10 ;
sol = pdepe(m,@pdepedx1,@pdepeic,@pdepebc,x,t);
u = sol(:,:,1);

figure
subplot 121
surf(x,t,u)
colorbar;
colormap jet
title('Numerical solution computed with 11 mesh points.',...
    'fontsi',20,'interpreter','latex')
xlabel('x','fontsi',30,'interpreter','latex')
ylabel('Time','fontsi',30,'interpreter','latex')
zlabel('c','fontsi',30,'interpreter','latex')