% Geomesh info
L = 0.4;                        % length of domain in x direction [m]       
tmax = 35;                      % end time [day]
nx = 41;                        % number of nodes in x direction
nt = 50401;                     % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);

I0 = 0.5 % Initial acid conc [M]


% Constant fractions
F = 96485;                        % Faraady constant [C/mol]
R = 8.314;                        % Gas constant [J/mol.K]


% Physical info
T = 25 + 273;                     % Temperature [K]
n = 0.64;                         % Porosity
tau = 1.25;                       % Tortuosity
z1 = 1;                           % Valancy constant (H+)
z2 = 0;                           % valancy constant (Hydrocarbon)
z3 = -1;                          % valancy constant (OH-)
sigma_surface = 2.74*10^7;        % Conductivity [S/m]



% Inputs which might change based on expriment
dEdx = 150;                       % Eletric field [V/m]
epsilon = 2.1;                    % dielectric constant [hydrocarbon]
epsilon_H = 1.2;                  % dielectric constant [Hydrogen]
epsilon_water = 80;               % dielectric constant [water]
epsilon_OH = 2.21;                % dielectric constant [OH]
zeta = -0.0027;                   % Zeta potential [V]
T = 25+273;                       % Absolute temperature [K]

mu_water = 0.001*24*3600;         % Water Viscosity [kg/(m.day)]
mu_oil = 510*24*3600;             % Oil viscosity [kg/(m.day)]
K = 0.02;


% H+ Properties
D_H = 3.5447*10^-9*24*3600;       % Mass advection [m^2/day]
D_star_H = D_H*coeff*(dt/dx^2);   % Dimensionless of Diffusion
D_star1_H = D_H*(dt/dx^2);
alpha_H = D_star_H/n;             % Diffusion Advection 
alpha1_H = D_star1_H/n;


% Hydrocarbon Properties
D = 2.063*10^-9*24*3600;          % Mass advection [m^2/day]
D_star = D*coeff*(dt/dx^2);       % Dimensionless of Diffusion
D_star1 = D*(dt/dx^2);
alpha = D_star/n;                 % Diffusion Advection 
alpha1 = D_star1/n;


% OH- Properties
D_OH = 0.450*10^-8*24*3600;       % Mass advection [m^2/day]
D_star_OH = D_OH*coeff*(dt/dx^2); % Dimensionless of Diffusion with K influence
D_star1_OH = D_OH*(dt/dx^2);      % Dimensionless of Diffusion without influence
alpha_OH = D_star_OH/n;           % Diffusion Advection 
alpha1_OH = D_star1_OH/n;
u_eo_OH = (epsilon_OH*zeta)/mu_water
%--------------------------------------------------------------------

% Electroosmotic Values
u_star = ((epsilon*zeta)/mu_oil)*(dEdx/(tau^2));                 % Hydrocarbon
u_star_H = ((epsilon_H*zeta)/mu_water)*(dEdx/(tau^2));           % H+
u_star_water = ((epsilon_water*zeta)/mu_water)*(dEdx/(tau^2));   % Water
u_star_OH = ((epsilon_OH*zeta)/mu_water)*(dEdx/(tau^2));         % OH-

%--------------------------------------------------------------------

% Electromigration values for Hydrocarbon
v = (D/(R*T));                            % mobility [Hydrocarbon]
u_e = (v*z2*F*dEdx)/(tau^2);              % electromigration [Hydrocarbon]
u_total = coeff*(u_e - u_star);           % Total mobility advection [Hydrocarbon]
u_total1 = (u_e - u_star);                % Total mobility advection without influence [Hydrocarbon]
beta = (u_total/n)*(dt/2*dx);             % Dimensionless mobility advection [Hydrocarbon]
beta1_H = (u_total1/n)*(dt/2*dx);         % Dimensionless mobility advection without K influence [Hydrocarbon]


% Electromigration values for H+
v_H = (D_H/(R*T));                        % mobility [H+]
u_e_H = (v_H*z1*F*dEdx)/(tau^2);          % electromigration [H+]
u_total_H = coeff*(u_e_H - u_star_H);     % Total mobility advection with K influence [H+]
u_total1_H = (u_e_H - u_star_H);          % Total mobility advection without influence [H+]
beta_H = (u_total_H/n)*(dt/2*dx);         % Dimensionless mobility advection with K influence [H+]
beta1_H = (u_total1_H/n)*(dt/2*dx);       % Dimensionless mobility advection without K influence [H+]


% Electromigration values for OH+
v_OH = (D_OH/(R*T));                      % mobility [OH-]
u_e_OH = (v_OH*z1*F*dEdx)/(tau^2);        % electromigration [OH-]
u_total_OH = coeff*(u_e_OH - u_star_OH);  % Total mobility advection [OH-]
u_total1_OH = (u_e_OH - u_star_OH);       % Total mobility advection without influence
beta_OH = (u_total_OH/n)*(dt/2*dx);       % Dimensionless mobility advection [OH-]
beta1_H = (u_total1_OH/n)*(dt/2*dx);      % Dimensionless mobility advection without K influence [H+]


% Flux
J0 = I0/sqrt(R1*alpha); 
J0 = J0/100;
% used for boundary condition



% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);

Sigma = zeros(nx,nt):
Sigma(:,1) = (F^2)*(z1^2)*(v_H)*I0;
Sigma(:,2) = (F^2)*(z1^2)*(v_H)*I0;

ix_H = zeros(nx,nt);


i = -Sigma*dEdx - 
% for both i and sigma we are in need of concetration so
% it will Calculated >>>>
