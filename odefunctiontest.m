clc,clear
m = 0;
x = linspace(0,0.4,20);
t = linspace(0,2,20);


sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);

u = sol;
surf(x,t,u)
title('Numerical solution computed with 50 mesh points.')
xlabel('Distance x')
ylabel('Time t')

figure
plot(x,u(end,:))
title('Solution at t = 2')
xlabel('Distance x')
ylabel('u(x,2)')

% local function 1
function [c,f,s] = pdex1pde(x,t,u,DuDx)

%Mian values
n = 0.64; %porosity
D0 = 1*10^-9; %refrence diffusivity
D = 1.2; %Diffusion (for i componet)
Tau = 1.25; %tortuosity
R_gas = 8.314; %Gas constant [J/mol.k-1]
Temp_cel = 15;
T = Temp_cel + 273; %Temp in Kelvin

%dimensesless values
Pe = 47;
Z = 0.049;
Beta = 967;

%dimenssial values
zeta = -2.7; %zeta potential
sigma_surface = 0.0013; %surface conductivity
Voltage = 25; %applied voltage
L = 0.4; %Length of porous medium
Ka = 1.75*10^-5; %Dissociation constant
NaCl_int = 0.003; % NaCl initial concetration
c_Na = 0.1; %purge solotion NaCl conc
c_cl = 0.1; %purge solotion NaCl conc
mu = 0.001; %Viscosity
epsilon = 7*10^-10; %Electrical permitivity

u_c = (1/Tau^2)*(zeta*Voltage);
u_s = n*u_c;
V = D/(R_gas*T); %Mobility
F = 9.64853321233100184*10^4; %Faraday Constant [C/mol]
z = 3; %charge number of componet
u_e = -V*z*F*Voltage*(1/Tau^2);
u_t = (u_e + u_c);
R = 0; %This is molar production which in this case it is zero
adsorption_coeff = 0.087;



c = 1;
f = D / Tau^2 * DuDx;
s = DuDx*1 + R;
end

% local function 2
function u0 = pdex1ic(x)
u0 = 0.5;
end

% local function 3
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;
end
