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

% Refrence t directions        [m]
t = (0:dt:(nt-1)*dt);
t_ref = repmat(t,nx,1);
t_up = t_ref/tmax;
t_step = t_up(1,2) - t_up(1,1);
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
K_H2O = 10^-8;                 % dissociation constant [(mol/m3)2]
K_b = 1.75*10^-6;              % dissociation constant [mol/m3]
mu_a = 0.001;                  % Solution viscosity [kg/m.s]
epsilon = 7*10^10;             % Electrical permittivity [F/m]
zeta = -0.0027;                % Zeta Potential [V]
zeta_0 = 2.6205e-23;           % Refrence Zeta Potential [V]

K_ads = 0.02;
coeff = 1/(1+K_ads);

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
x_step = 0.250;
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
u_c = u_x/((tau^2)*10^19);              % Convection velocity
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
u_t_A_up = (u_e_A_upp + (u_c_up));
u_t_Na_up = (u_e_Na_upp + u_c_up);
u_t_Cl_up = (u_e_Cl_upp + u_c_up);
u_t_H_up = (u_e_H_upp + u_c_up)*10^-5;
u_t_OH_up = (u_e_OH_upp + u_c_up);
u_t_C_up = (u_e_C_upp + u_c_up);

% Initial concentration        [mol/m3]
c_0 = 500; 
c_p = 200;
c_Na = c_p;
c_Cl = c_Na;

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

% concentration arrays [remapped]
G_HA_up = zeros(nx,nt);
G_A_up = zeros(nx,nt);
G_Na_up = zeros(nx,nt);
G_Cl_up = zeros(nx,nt);
G_H_up = zeros(nx,nt);
G_OH_up = zeros(nx,nt);
G_C_up = zeros(nx,nt);

% adsorbed concentration arrays [remapped]
G_HA_ads = zeros(nx,nt);
G_A_ads = zeros(nx,nt);
G_Na_ads = zeros(nx,nt);
G_Cl_ads = zeros(nx,nt);
G_H_ads = zeros(nx,nt);
G_OH_ads = zeros(nx,nt);
G_C_ads = zeros(nx,nt);

% Current array
i_z = zeros(nx, nt);

% Sigma array
sigma_bar = zeros(nx,nt);
sum_HA = zeros(nx,nt);
sum_A = zeros(nx,nt);
sum_Na = zeros(nx,nt);
sum_Cl = zeros(nx,nt);
sum_H = zeros(nx,nt);
sum_OH = zeros(nx,nt);
sum_C = zeros(nx,nt);
sum_total = zeros(nx,nt);

C_HA = zeros(nx,nt);
C_A = zeros(nx,nt);
C_Na = zeros(nx,nt);
C_Cl = zeros(nx,nt);
C_H = zeros(nx,nt);
C_OH = zeros(nx,nt);
C_C = zeros(nx,nt);
C_total = zeros(nx,nt);

% SpeciesConstants
K_H2O_m = ones(nx,nt);
K_H2O_m = K_H2O_m*K_H2O;

K_a_m = ones(nx,nt);
K_a_m = K_a_m*K_a;

K_b_m = ones(nx,nt);
K_b_m = K_b_m*K_b;

% Rate arrays
R_HA = zeros(nx,nt);
R_A = zeros(nx,nt);
R_Na = zeros(nx,nt);
R_Cl = zeros(nx,nt);
R_H = zeros(nx,nt);
R_OH = zeros(nx,nt);
R_C = zeros(nx,nt);

%Adsorbed Rate arrays
R_HA_ads = zeros(nx,nt);
R_A_ads = zeros(nx,nt);
R_Na_ads = zeros(nx,nt);
R_Cl_ads = zeros(nx,nt);
R_H_ads = zeros(nx,nt);
R_OH_ads = zeros(nx,nt);
R_C_ads = zeros(nx,nt);

% --- Set IC and BC
G_HA(:,1) = c_0;
G_A(:,1) = c_0;
G_Na(:,1) = c_Na;
G_Cl(:,1) = c_Cl;
G_H(:,1) = c_0;
G_OH(:,1) = c_0;
G_C(:,1) = c_0;

G_HA_up(:,1) = G_HA(:,1)/c_0;
G_A_up(:,1) = G_A(:,1)/c_0;
G_Na_up(:,1) = G_Na(:,1)/c_Na;
G_Cl_up(:,1) = G_Cl(:,1)/c_Cl;
G_H_up(:,1) = G_H(:,1)/c_0;
G_OH_up(:,1) = G_OH(:,1)/c_0;
G_C_up(:,1) = G_C(:,1)/c_0;

G_HA_ads(:,1) = K_ads*c_0;
G_A_ads(:,1) = K_ads*c_0;
G_Na_ads(:,1) = K_ads*c_Na;
G_Cl_ads(:,1) = K_ads*c_Cl;
G_H_ads(:,1) = K_ads*c_0;
G_OH_ads(:,1) = K_ads*c_0;
G_C_ads(:,1) = K_ads*c_0;

% Fixing alpha advection
alpha_HA = (D_HA/(tau^2))*10^5;
alpha_A = (D_A/(tau^2))*10^5;
alpha_Na = (D_Na/(tau^2))*10^5;
alpha_Cl = (D_Cl/(tau^2))*10^5;
alpha_H = (D_H/(tau^2))*10^5;
alpha_OH = (D_OH/(tau^2))*10^5;
alpha_C = (D_C/(tau^2))*10^5;

for m=1:nt-1
    % --- Set IC and BC
    G_HA_up(1,m) = u_t_HA_up(1,m);
    G_HA_up(end,m) = 0;

    G_A_up(1,m) = u_t_A_up(1,m);
    G_A_up(end,m) = 0;

    G_Na_up(1,m) = u_t_Na_up(1,m);
    G_Na_up(end,m) = 0;

    G_Cl_up(1,m) = u_t_Cl_up(1,m);
    G_Cl_up(end,m) = 0;

    G_H_up(1,m) = u_t_H_up(1,m);
    G_H_up(end,m) = 0;        

    G_OH_up(1,m) = u_t_OH_up(1,m);
    G_OH_up(end,m) = 0;

    G_C_up(1,m) = u_t_C_up(1,m);
    G_C_up(end,m) = 0;
    % ---
    G_HA_ads(1,m) = u_t_H(1,m)*c_0;
    G_HA_ads(end,m) = 0;
    
    G_A_ads(1,m) = u_t_H(1,m)*c_0;
    G_A_ads(end,m) = 0;  
    
    G_Na_ads(1,m) = u_t_H(1,m)*c_Na;
    G_Na_ads(end,m) = 0;  
    
    G_Cl_ads(1,m) = u_t_H(1,m)*c_Cl;
    G_Cl_ads(end,m) = 0;  

    G_H_ads(1,m) = u_t_H(1,m)*c_0;
    G_H_ads(end,m) = 0;  

    G_OH_ads(1,m) = u_t_H(1,m)*c_0;
    G_OH_ads(end,m) = 0;

    G_C_ads(1,m) = u_t_H(1,m)*c_0;
    G_C_ads(end,m) = 0;  

    for i=2:nx-1
        % Adsorbed concentration is used for rate of reaction in each term  
        G_H_ads(i,m+1) = G_H_ads(i,m) + coeff*(alpha_H*((G_H_ads(i+1,m) -2*G_H_ads(i,m) + G_H_ads(i-1,m))*((dt/n)/(dx^2))) - (((u_t_H(i+1,m) - u_t_H(i,m))*((dt/n)/dx)) * ((G_H_ads(i+1,m) - G_H_ads(i,m))*((dt/n)/dx))));
        G_A_ads(i,m+1) = G_A_ads(i,m) + coeff*(alpha_A*((G_A_ads(i+1,m) -2*G_A_ads(i,m) + G_A_ads(i-1,m))*((dt/n)/(dx^2))) - (((u_t_A(i+1,m) - u_t_A(i,m))*((dt/n)/dx)) * ((G_A_ads(i+1,m) - G_A_ads(i,m))*((dt/n)/dx))));
        G_C_ads(i,m+1) = G_C_ads(i,m) + coeff*(alpha_C*((G_C_ads(i+1,m) -2*G_C_ads(i,m) + G_C_ads(i-1,m))*((dt/n)/(dx^2))) - (((u_t_C(i+1,m) - u_t_C(i,m))*((dt/n)/dx)) * ((G_C_ads(i+1,m) - G_C_ads(i,m))*((dt/n)/dx))));
        G_HA_ads(i,m+1) = G_HA_ads(i,m) + coeff*(alpha_HA*((G_HA_ads(i+1,m) -2*G_HA_ads(i,m) + G_HA_ads(i-1,m))*((dt/n)/(dx^2))) - (((u_t_HA(i+1,m) - u_t_HA(i,m))*((dt/n)/dx)) * ((G_HA_ads(i+1,m) - G_HA_ads(i,m))*((dt/n)/dx))));
        G_Na_ads(i,m+1) = G_Na_ads(i,m) + coeff*(alpha_Na*((G_Na_ads(i+1,m) -2*G_Na_ads(i,m) + G_Na_ads(i-1,m))*((dt/n)/(dx^2))) - (((u_t_Na(i+1,m) - u_t_Na(i,m))*((dt/n)/dx)) * ((G_Na_ads(i+1,m) - G_Na_ads(i,m))*((dt/n)/dx))));
        G_Cl_ads(i,m+1) = G_Cl_ads(i,m) + coeff*(alpha_Cl*((G_Cl_ads(i+1,m) -2*G_Cl_ads(i,m) + G_Cl_ads(i-1,m))*((dt/n)/(dx^2))) - (((u_t_Cl(i+1,m) - u_t_Cl(i,m))*((dt/n)/dx)) * ((G_Cl_ads(i+1,m) - G_Cl_ads(i,m))*((dt/n)/dx))));
        G_OH_ads(i,m+1) = G_OH_ads(i,m) + coeff*(alpha_OH*((G_OH_ads(i+1,m) -2*G_OH_ads(i,m) + G_OH_ads(i-1,m))*((dt/n)/(dx^2))) - (((u_t_OH(i+1,m) - u_t_OH(i,m))*((dt/n)/dx)) * ((G_OH_ads(i+1,m) - G_OH_ads(i,m))*((dt/n)/dx))));
        
        R_H_ads(i,m+1) = G_H_ads(i,m+1);
        R_A_ads(i,m+1) = G_A_ads(i,m+1);
        R_C_ads(i,m+1) = G_C_ads(i,m+1);
        R_HA_ads(i,m+1) = G_HA_ads(i,m+1);
        R_Na_ads(i,m+1) = G_Na_ads(i,m+1);
        R_Cl_ads(i,m+1) = G_Cl_ads(i,m+1);
        R_OH_ads(i,m+1) = G_OH_ads(i,m+1);

        % the 10^5 is a fixing factor in H and can be different in other Species
        R_H(i,m+1) = R_H_ads(i,m+1)/((c_0*k_0)*10000);
        R_A(i,m+1) = R_A_ads(i,m+1)/((c_0*k_0)*10000);
        R_C(i,m+1) = R_C_ads(i,m+1)/((c_0*k_0)*10000);
        R_HA(i,m+1) = R_HA_ads(i,m+1)/((c_0*k_0)*10000);
        R_Na(i,m+1) = R_Na_ads(i,m+1)/((c_0*k_0)*10000);
        R_Cl(i,m+1) = R_Cl_ads(i,m+1)/((c_0*k_0)*10000);
        R_OH(i,m+1) = R_OH_ads(i,m+1)/((c_0*k_0)*10000);
        
        G_HA_up(i,m+1) = G_HA_up(i,m) + (D_HA_upp/(Peclet_calculated))*((G_HA_up(i+1,m) -2*G_HA_up(i,m) + G_HA_up(i-1,m))*((t_step/n)/(x_step^2))) - (((u_t_HA_up(i+1,m) - u_t_HA_up(i,m))*((t_step/n)/x_step)) * ((G_HA_up(i+1,m) - G_HA_up(i,m))*((t_step/n)/x_step)));
        G_Na_up(i,m+1) = G_Na_up(i,m) + (D_Na_upp/(Peclet_calculated))*((G_Na_up(i+1,m) -2*G_Na_up(i,m) + G_Na_up(i-1,m))*((t_step/n)/(x_step^2))) - (((u_t_Na_up(i+1,m) - u_t_Na_up(i,m))*((t_step/n)/x_step)) * ((G_Na_up(i+1,m) - G_Na_up(i,m))*((t_step/n)/x_step)));
        G_Cl_up(i,m+1) = G_Cl_up(i,m) + (D_Cl_upp/(Peclet_calculated))*((G_Cl_up(i+1,m) -2*G_Cl_up(i,m) + G_Cl_up(i-1,m))*((t_step/n)/(x_step^2))) - (((u_t_Cl_up(i+1,m) - u_t_Cl_up(i,m))*((t_step/n)/x_step)) * ((G_Cl_up(i+1,m) - G_Cl_up(i,m))*((t_step/n)/x_step)));
        G_OH_up(i,m+1) = G_OH_up(i,m) + (D_OH_upp/(Peclet_calculated))*((G_OH_up(i+1,m) -2*G_OH_up(i,m) + G_OH_up(i-1,m))*((t_step/n)/(x_step^2))) - (((u_t_OH_up(i+1,m) - u_t_OH_up(i,m))*((t_step/n)/x_step)) * ((G_OH_up(i+1,m) - G_OH_up(i,m))*((t_step/n)/x_step)));
        G_A_up(i,m+1) = G_A_up(i,m) + (D_A_upp/(Peclet_calculated))*((G_A_up(i+1,m) -2*G_A_up(i,m) + G_A_up(i-1,m))*((t_step/n)/(x_step^2))) - (((u_t_A_up(i+1,m) - u_t_A_up(i,m))*((t_step/n)/x_step)) * ((G_A_up(i+1,m) - G_A_up(i,m))*((t_step/n)/x_step)));
        G_H_up(i,m+1) = G_H_up(i,m) + (D_H_upp/(Peclet_calculated))*((G_H_up(i+1,m) -2*G_H_up(i,m) + G_H_up(i-1,m))*((t_step/n)/(x_step^2))) - (((u_t_H_up(i+1,m) - u_t_H_up(i,m))*((t_step/n)/x_step)) * ((G_H_up(i+1,m) - G_H_up(i,m))*((t_step/n)/x_step))) + (1/alpha(i,m+1))*R_H(i,m+1);        
        G_C_up(i,m+1) = G_C_up(i,m) + (D_C_upp/(Peclet_calculated))*((G_C_up(i+1,m) -2*G_C_up(i,m) + G_C_up(i-1,m))*((t_step/n)/(x_step^2))) - (((u_t_C_up(i+1,m) - u_t_C_up(i,m))*((t_step/n)/x_step)) * ((G_C_up(i+1,m) - G_C_up(i,m))*((t_step/n)/x_step)));
         
        
        sum_H(i,m) = z_H*D_H_upp*G_H_up(i,m);
        sum_A(i,m) = z_A*D_H_upp*G_A_up(i,m);
        sum_C(i,m) = z_C*D_H_upp*G_C_up(i,m);
        sum_HA(i,m) = z_HA*D_H_upp*G_HA_up(i,m);
        sum_Na(i,m) = z_Na*D_H_upp*G_Na_up(i,m);
        sum_Cl(i,m) = z_Cl*D_H_upp*G_Cl_up(i,m);
        sum_OH(i,m) = z_OH*D_H_upp*G_OH_up(i,m);

        % Sigma calculations
        sum_total(i,m) = sum_H(i,m) + sum_A(i,m) + sum_C(i,m) + sum_HA(i,m) + sum_Na(i,m) + sum_Cl(i,m) + sum_OH(i,m);
        sigma_bar(i,m) = (1/(tau^2))*sum_total(i,m);

        C_H (i,m) = (G_H_up(i+1,m) - G_H_up(i,m))*((t_step/n)/x_step);
        C_A (i,m) = (G_A_up(i+1,m) - G_A_up(i,m))*((t_step/n)/x_step);
        C_C (i,m) = (G_C_up(i+1,m) - G_C_up(i,m))*((t_step/n)/x_step);
        C_HA (i,m) = (G_HA_up(i+1,m) - G_HA_up(i,m))*((t_step/n)/x_step);
        C_Na (i,m) = (G_Na_up(i+1,m) - G_Na_up(i,m))*((t_step/n)/x_step);
        C_Cl (i,m) = (G_Cl_up(i+1,m) - G_Cl_up(i,m))*((t_step/n)/x_step);
        C_OH (i,m) = (G_OH_up(i+1,m) - G_OH_up(i,m))*((t_step/n)/x_step);
        C_total(i,m) = (D_H_upp*z_H*C_H(i,m) + D_A_upp*z_A*C_A(i,m) + D_C_upp*z_C*C_C(i,m) + D_HA_upp*z_H*C_HA(i,m) + D_Na_upp*z_Na*C_Na(i,m) + D_Cl_upp*z_Cl*C_Cl(i,m) + D_OH_upp*z_OH*C_OH(i,m));
        
        i_z(i,m) = sigma_bar(i,m)*dphidx(i,m) - (1/Beta_calculated)*(C_total(i,m));
        
        R_H(1,m) = (i_z(2,m)/F)*1000;
        R_OH(1,m) = (i_z(2,m)/F)*1000;

        for ii=1:nx
            if m == 1
                R_H(ii,m) = (i_z(2,m)/F)*(Beta/(Peclet*Z)) + 0.1;
            end
        end
      
    end
end





pH = log10(G_H);

pH_scale = linspace(1,40,41);
xl = [0,5,10,15,20,25,30,35];
yl = [10000,7900,7100,6000,5700,5500,5400,5100];

figure(1);

hold on;


%plot(t_array,G_HA_up(10,:),'-','DisplayName', 'HA');
%plot(t_array,G_A_up(10,:),'-','DisplayName', 'A-');
%plot(t_array,G_Na_up(10,:),'-','DisplayName', 'Na+');
%plot(t_array,G_Cl_up(10,:),'-','DisplayName', 'Cl'); $
plot(t_array,G_H_up(10,:),'-','DisplayName', 'H+ ');
%plot(t_array,R_H(10,:),'-','DisplayName', 'H+++ ');
%plot(t_array,G_OH_up(10,:),'-','DisplayName', 'OH- ');
%plot(t_array,G_C_up(10,:),'-','DisplayName', 'Hydrocarbon');

legend();
