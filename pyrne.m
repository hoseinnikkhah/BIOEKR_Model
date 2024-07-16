% Geomesh info
L = 0.4;                       % length of domain in x direction [m]
L_cm = L*100;                  % length of domain in x direction [cm]

% Nodes
tmax = 35;                     % end time [day]
nx = 41;                       % number of nodes in x direction
nt = 50401;                    % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);

% Refrence x directions        [m]
x = (10^-5:dx:(nx)*dx);
x = transpose(x);
x_ref = repmat(x,1,nt);

% Physical info 
V = 25;                        % Voltage [V]

E_field = ones(nx,nt);         % Voltage [V]
M = linspace(V,0,nx);
for timestep = 1:nt
    E_field(:,timestep) = M;
end

dVdx = V/L;                    % Voltage gradient [V/m]
dVdx_cm = V/L_cm;              % Voltage gradient [V/cm]
T = 37 + 273;                  % Tempature [K]
F = 96485;                     % Faraady constant [C/mol]
R = 8.314;                     % Gas constant [J/mol.K]
D0 = 10^-9;                    % Reference diffusivity [m2/s]

% Soil info
n = 0.64;                      % Porosity   [Dimentionless]
tau = 1.25;                    % Tortuosity [Dimentionless]

% Acetic acid
sigma_surface = 0.0013;        % Surface conductivit [mhos/m]
K_a = 1.75*10^-6;              % dissociation constant [mol/m3]
mu_a = 0.001;                  % Solution viscosity [kg/m.s]
epsilon = 7*10^10;             % Electrical permittivity [F/m]
zeta = -0.0027;                % Zeta Potential [V]
zeta_0 = 2.6205e-23;           % Refrence Zeta Potential [V]

coeff = 1/(1+K_a);
% Dimentionless     [Dimentionless]
Peclet = 47;
Z = 0.049;
Beta = 967;
k_0 = K_a;

% Valecy            [Dimentionless]
z_HA = 0;
z_A = -1;
z_Na = 1;
z_Cl = -1;
z_H = 1;
z_OH = -1;
z_C = 0;
% In an array
z_i = [z_HA, z_A, z_Na, z_Cl, z_H, z_OH,z_C];

% Species diffusivities (remapped)  [Dimentionless]
D_HA_upp = 1.2;          % Acetic Acid
D_A_upp = 1.2;           % Acid Agent
D_Na_upp = 1.34;         % Na+
D_Cl_upp = 2.05;         % Cl-
D_H_upp = 9.35;          % H+
D_OH_upp = 2.00;         % OH-
D_C_upp = 2.00;          % Carbon

D_i_upp = [D_HA_upp, D_A_upp, D_Na_upp, D_Cl_upp, D_H_upp, D_OH_upp, D_C_upp];

% Species diffusivities (Normal)    [m2/s]
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

% Electroelectromigration velocity (Normal)     [m/s]
u_e_HA = -v_HA*z_HA*F*E_field*(1/tau^2);
u_e_A = -v_A*z_A*F*E_field*(1/tau^2);
u_e_Na = -v_Na*z_Na*F*E_field*(1/tau^2);
u_e_Cl = -v_Cl*z_Cl*F*E_field*(1/tau^2);
u_e_H = -v_H*z_H*F*E_field*(1/tau^2);
u_e_OH = -v_OH*z_OH*F*E_field*(1/tau^2);
u_e_C = -v_C*z_C*F*E_field*(1/tau^2);

% Refrence velocity                             [m/s]
u_0 = (1/tau^2)*((epsilon*zeta)/mu_a)*E_field;

% Dimensionless x directions
x_up = x_ref/L;

% Dimensionless Voltage     [Dimentionless]
phi_bar = E_field/V;
dphidx = phi_bar./x_up;

% Dimentionless Calculated     [Dimentionless]
Peclet_calculated = (epsilon*zeta_0*V)/(mu_a*D0);
Beta_calculated = (F*V)/(R*T);
Z_calculated = (R*T*epsilon*zeta_0)/(D0*F*mu_a);

% Electroelectromigration velocity (remapped)        [Dimentionless]
u_e_HA_upp = (-1/Z_calculated)*D_HA_upp*z_HA*dphidx;
u_e_A_upp = (-1/Z_calculated)*D_A_upp*z_A*dphidx;
u_e_Na_upp = (-1/Z_calculated)*D_Na_upp*z_Na*dphidx;
u_e_Cl_upp = (-1/Z_calculated)*D_Cl_upp*z_Cl*dphidx;
u_e_H_upp = (-1/Z_calculated)*D_H_upp*z_H*dphidx;
u_e_OH_upp = (-1/Z_calculated)*D_OH_upp*z_OH*dphidx;
u_e_C_upp = (-1/Z_calculated)*D_C_upp*z_C*dphidx;


% Convection velocity          [m/s]
u_x = (epsilon/mu_a)*(zeta*E_field);
u_c = u_x/((tau^2)*10^14);              % Convection velocity
u_c_up = -((zeta/zeta_0)*dphidx)*10^-18;

% Toatal velocity term (Normal)
u_t_HA = (u_e_HA + u_c);
u_t_A = (u_e_A + u_c);
u_t_Na = (u_e_Na + u_c);
u_t_Cl = (u_e_Cl + u_c);
u_t_H = (u_e_H + u_c);
u_t_OH = (u_e_OH + u_c);
u_t_C = (u_e_C + u_c);

% Toatal velocity term (Rmapped)
u_t_HA_up = (u_e_HA_upp + u_c_up);
u_t_A_up = (u_e_A_upp + u_c_up);
u_t_Na_up = (u_e_Na_upp + u_c_up);
u_t_Cl_up = (u_e_Cl_upp + u_c_up);
u_t_H_up = (u_e_H_upp + u_c_up);
u_t_OH_up = (u_e_OH_upp + u_c_up);
u_t_C_up = (u_e_C_upp + u_c_up);

% Initial concentration        [mol/m3]
c_0 = 500; 

% Sigma refrence               [S/m]
sigma_0 = ((F^2)*D0*c_0*(tau^2))/(R*T);

% Alpha advection
alpha = u_0/(c_0*k_0);


% --- Create arrays to save data for export
x_array = linspace(0,L,nx);
t_array = linspace(0,tmax,nt);

% Flux arrays
J_HA = zeros(nx,nt);
J_A = zeros(nx,nt);
J_Na = zeros(nx,nt);
J_Cl = zeros(nx,nt);
J_H = zeros(nx,nt);
J_OH = zeros(nx,nt);
J_C = zeros(nx,nt);

% concentration arrays
G_HA = zeros(nx,nt);
G_A = zeros(nx,nt);
G_Na = zeros(nx,nt);
G_Cl = zeros(nx,nt);
G_H = zeros(nx,nt);
G_OH = zeros(nx,nt);
G_C = zeros(nx,nt);

% Sigma arrays
Sigma = zeros(nx,nt);
Sigma_ref = ones(nx,nt);
sigma_ref = Sigma_ref*sigma_surface;

% Current array
i_z = zeros(nx, nt);

% species sigma
s_H = zeros(nx,nt);
s_C = zeros(nx,nt);
s_OH = zeros(nx,nt);

% SpeciesConstants
K_H2O = zeros(nx,nt);
K_a = zeros(nx,nt);
K_b = zeros(nx,nt);

% Rate arrays
R_HA = zeros(nx,nt);
R_A = zeros(nx,nt);
R_Na = zeros(nx,nt);
R_Cl = zeros(nx,nt);
R_H = zeros(nx,nt);
R_OH = zeros(nx,nt);
R_C = zeros(nx,nt);

% --- Set IC and BC
%G_HA(:,1) = c_0;
% J_H(1,:) = u_t_H(1,:)*c_0;

%for m=1:nt-1

%    for i=2:nx-1
%        G_HA(1,m) = (u_t_H(1,m)*c_0) - (D_HA/(tau^2))*((G_HA(i,m) - G_HA(i-1,m))*(dt/(n*dx)));        %--- z = 0 boundary
%        G_HA(end,m) = (u_t_H(1,m)*c_0) - (D_HA/(tau^2))*((G_HA(i,m) - G_HA(i-1,m))*(dt/(n*dx)));      %--- z = L boundary
%        G_HA(i,m+1) = G_HA(i,m) + coeff*((D_HA/(tau^2))*((G_HA(i+1,m) -2*G_HA(i,m) + G_HA(i-1,m))*(dt/(n*(dx^2)))) - (((u_t_HA(i+1,m) - u_t_HA(i-1,m))*(dt/(n*2*dx))) * ((G_HA(i+1,m) - G_HA(i-1,m))*(dt/(n*2*dx)))));
        %G_H(i,m+1) = G_H(i,m) + coeff*((D_H/(tau^2))*((G_H(i+1,m) -2*G_H(i,m) + G_H(i-1,m))*(dt/(n*(dx^2)))) - ((u_t_H(i+1,m) - u_t_H(i-1,m))*(G_H(i+1,m) - G_H(i-1,m)))*(dt/(n*2*dx)));
        
%    end
%end




G_HA_up = zeros(nx,nt);

for m=1:nt-1
    G_HA_up(1,m) = 1;
    G_HA_up(end,m) = 0;
    for i=2:nx-1
        G_HA_up(i,m+1) = G_HA_up(i,m) + (D_HA_upp/(Peclet_calculated))*((G_HA_up(i+1,m) -2*G_HA_up(i,m) + G_HA_up(i-1,m))*(dt/((x_up(i,m))^2))) - (((u_t_HA(i+1,m) - u_t_HA(i-1,m))*(dt/(n*2*dx))) * ((G_HA(i+1,m) - G_HA(i-1,m))*(dt/(n*2*dx))));
    end
end

pH = log10(G_H);

pH_scale = linspace(1,40,41);
xl = [0,5,10,15,20,25,30,35];
yl = [10000,7900,7100,6000,5700,5500,5400,5100];
figure;

hold on;

plot(t_array,G_HA_up(10,:),'-','DisplayName', 'Hydrocarbon (EKR)');

