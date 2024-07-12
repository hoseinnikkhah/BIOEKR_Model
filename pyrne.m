% Geomesh info
L = 0.4;                       % length of domain in x direction [m]
L_cm = L*100;                  % length of domain in x direction [cm]

% Nodes
tmax = 35;                     % end time [day]
nx = 41;                       % number of nodes in x direction
nt = 50401;                    % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);


% Physical info 
V = 25;                        % Voltage [V]
dVdx = V/L;                    % Voltage gradient [V/m]
dVdx_cm = V/L_cm;              % Voltage gradient [V/cm]
T = 37 + 273;                  % Tempature [K]
F = 96485;                     % Faraady constant [C/mol]
R = 8.314;                     % Gas constant [J/mol.K]
D0 = 10^-9;                    % Reference diffusivity [m2/s]

% Soil info
n = 0.64;                      % Porosity
tau = 1.25;                    % Tortuosity

% Acetic acid
sigma_surface = 0.0013;        % Surface conductivit [mhos/m]
K_a = 1.75/100000;             % dissociation constant [mol/m3]
mu_a = 0.001;                  % Solution viscosity [kg/m.s]
epsilon = 7*10^10;             % Electrical permittivity [F/m]
zeta = -0.0027;                % Zeta Potential [V]

% Dimentionless 
Peclet = 47;
Z = 0.049;
Beta = 967;

% Valecy
z_HA = 0;
z_A = -1;
z_Na = 1;
z_Cl = -1;
z_H = 1;
z_OH = -1;
z_C = 0;
% In an array
z_i = [z_HA, z_A, z_Na, z_Cl, z_H, z_OH,z_C];

% Species diffusivities (remapped)
D_HA_upp = 1.2;          % Acetic Acid
D_A_upp = 1.2;           % Acid Agent
D_Na_upp = 1.34;         % Na+
D_Cl_upp = 2.05;         % Cl-
D_H_upp = 9.35;          % H+
D_OH_upp = 2.00;         % OH-
D_C_upp = 2.00;          % Carbon

D_i_upp = [D_HA_upp, D_A_upp, D_Na_upp, D_Cl_upp, D_H_upp, D_OH_upp, D_C_upp];

% Species diffusivities (Normal)
D_HA = D_HA_upp*D0;      % Acetic Acid
D_A = D_A_upp*D0;        % Acid Agent
D_Na = D_Na_upp*D0;      % Na+
D_Cl = D_Cl_upp*D0;      % Cl-
D_H = D_H_upp*D0;        % H+
D_OH = D_OH_upp*D0;      % OH-
D_C = D_C_upp*D0;        % Carbon

D_i = [D_HA, D_A, D_Na, D_Cl, D_H, D_OH, D_C];

% Mobility (Normal)      % [s·mol/kg]
v_HA = D_HA/(R*T);
v_A = D_A/(R*T);
v_Na = D_Na/(R*T);
v_Cl = D_Cl/(R*T);
v_H = D_H/(R*T);
v_OH = D_OH/(R*T);
v_C = D_C/(R*T);

v_i = [v_HA, v_A, v_Na, v_Cl, v_H, v_OH, v_C];

% Mobility (remapped)      % [s·mol/kg]
v_HA_upp = D_HA_upp/(R*T);
v_A_upp = D_A_upp/(R*T);
v_Na_upp = D_Na_upp/(R*T);
v_Cl_upp = D_Cl_upp/(R*T);
v_H_upp = D_H_upp/(R*T);
v_OH_upp = D_OH_upp/(R*T);
v_C_upp = D_C_upp/(R*T);

v_i_upp = [v_HA_upp, v_A_upp, v_Na_upp, v_Cl_upp, v_H_upp, v_OH_upp, v_C_upp];




