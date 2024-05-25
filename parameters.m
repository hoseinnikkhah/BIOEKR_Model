% Global m file for parameters need across domain

% Geometry parameters
L = 40;                           % length of domain in x direction [cm]
t_obs = 35;                       % total time [days]
tmax = t_obs*(24*60);             % total time [min]
nx = 41;                          % number of nodes in x direction
nt = tmax +1;                     % number of time steps
dx = L/(nx-1);                    % position step size [cm]
dt = tmax/(nt-1);                 % time step size [min]

% Global constants
R = 8.314;                        % Gas constant [J/mol.K]
F = 96485.33212;                  % Faraday constant [C/mol]

% Physical parameters (are given)
D0 = 3.5447*10^-9;                % Diffusion (Crud on in clay) [m^2/s] 
D1 = 9.3100*10^-8;                % Diffusion (H+ on in clay)
T = 25+273;                       % Absolute temperature [K]
mu_oil = 51;                      % Oil Viscosity [kg/(m.s)]
mu_water = 1;                     % Water Viscosity [kg/(m.s)]
dEdx = 1.5;                       % Voltage gradient [V/cm]
zeta = -0.0027;                   % Zeta potential [V]
c0 = 12;                          % Crude oil initial concentration [mg/kg]
c1 = 2;                           % Water initial concentration [mol/L]
u_e = 0;                          % Electromigration mobility [m^2/V.s] 
e_0 = 8.854*10^-12;               % permittivity of free space [F/m]
e_r = 7.5;                        % relative permittivity of clay [F/m]
e_oil = 2.3;                      % Crude oil permittivity [F/m]
e_clay = e_0*e_r;                 % Clay permittivity [F/m]
e_w = 5;                          % Water permittivity [F/m]
n=0.64;                           % Porosity
tau = 0.44;                       % Tortuosity
z = 1;                            % Valency
rho = 1620;                       % Clay Density [kg/m^3]
ol = 3.7*10^-9;                   % Conductivity for Ni in kaolinite [S/m]
% Calculated constants
k_eo1 = (e_oil*zeta)*n/(mu_oil);                 % Electroosmotic mobility Alshawabkeh [m^2/V.s]
k_eo2 = -(e_oil*zeta)*n/(mu_oil*(tau^2));        % Electroosmotic mobility Vane [m^2/V.s]
k_eo3 = (e_oil*zeta)*dEdx/(mu_oil);              % Electroosmotic mobility Shapiro [m^2/V.s]

Debye = sqrt((e_oil*T)/(2*(F^2)*(z^2)*c0));  % Debye length
J0 = -D0*c0 - (c0*(u_e + k_eo3)*dEdx);       % initial crude oil flux
Jw = -D1*c1 - (c1*(u_e + k_eo3)*dEdx);       % initial crude oil flux
V = L*dEdx;                                  % Applied Voltage [V]

R1 = - 0.5;
% --- Create arrays to save data for export

% --- Set IC and BC


