clc
%------------------------------------------------------------------------%
% global inputs
global nnd nel nne nodof eldof n %(active DOF)

%------------------------------------------------------------------------%
%these are parameters needed but may not be in equtations as a term
double_layer = 10; % double layer thickness in nm
a = 0.01; % porus radious
% in soils it is between 0.1 to 0.01 mm, the more dense the smaller a is

%------------------------------------------------------------------------%
%Dimensional Parameters for Acetic Acid Simulations
initial_conc_case1 = 0.5; % initial concentration of acetic acid in case 1
initial_conc_case2 = 0.1; % initial concentration of acetic acid in case 2
zeta = -2.7; % Zeta potential
F = 96485.33212; % C/mol
voltage = 25; % Applied voltage, in most cases this is constant
length = 0.4; % Porouse length which should be later in finite length
k_a = 1.75*10^-5; % Dissociation constant for acid
k_H2O = 10^-8; % Dissociation constant for water
initial_conc_NaCl = 0.003; % Initial NaCl concentration
NaCl_purge = 0.10; % Purge solution NaCl conentration (M)
mu = 0.001; % solution viscosity
epsilon = 7*10^-10; % Electrical permittivity
D0 = 10^-9; % Reference diffusivit
R = 8.314; % J/mol.k
T = 273 + 25;
% epsilon = [6.3 3.2 2.7 78.54*10^-22 1.000264 80]; ***********
z_i = [0 -1 +1 -1 +1 -1 +1 -1];
% for elements ch3cooh,ch3coo,Na,Cl,H and OH in order

%------------------------------------------------------------------------%
D_i = [1.2 1.2 1.34 2.05 9.35 2.00]; %m total is 6
% for elements ch3cooh,ch3coo,Na,Cl,H and OH in order they are seprated below
%check https://www.aqion.de/site/diffusion-coefficients for more info all
%the numbers must be devided to 10^-9
D_ch3cooh = 1.2;
D_ch3coo = 1.2;
D_Na = 1.34;
D_Cl = 2.05;
D_H = 9.35;
D_OH = 5.27;
v_i = D_i/T;
%------------------------------------------------------------------------%
%Sigma is surface and solution combined
sigma_surface = 0.0013; % Surface conductivity
% sigma can not be calculated like v_i, it is related to c_i and c_i is in
%one of pde problems so it should have its own function on a seprated .m
%file

%------------------------------------------------------------------------%


%------------------------------------------------------------------------%
%Nondimensional Parameters for Acetic Acid Simulations
np = 0.64; % Porosity
tau = 1.25; % Tortuosity
Pe = 47; % Pecklet number
Z = 0.049;
beta = 967;
Z1 = (R*T*epsilon*zeta)/(D0*F*mu);
Z1 = (-1*Z1)/1000;
pe = (epsilon*zeta*voltage)/(mu*D0);
pe = (-1*pe)/1000;
BETA = pe/Z1;
%------------------------------------------------------------------------%
%velocity terms
u_x = 7.5*10^-8; % average velocity in x direction
u_c = u_x/tau; % convenction velocity
u_s = np*u_c; % flow rate per unit area of the porous medium
%------------------------------------------------------------------------%
%here are constants that are valid all over length


