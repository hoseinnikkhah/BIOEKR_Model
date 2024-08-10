% Guide Book
% ---
% Each section is devided, and has 3 distinct states
% $$$ == fully sophiscated Check
% ### == not a sophiscated check
% !!! == unchecked
% These comments will be held on top and end of each section
% ---

% $$$
% Geomesh info
L = 0.4;                       % length of domain in x direction [m]

% $$$
% Nodes
tmax = 35;                     % end time [day]
nx = 41;                       % number of nodes in x direction
nt = 50401;                    % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);

% $$$
% Refrence x directions        [m]      $$$
x = (10^-5:dx:(nx)*dx);
x = transpose(x);
x_ref = repmat(x,1,nt);

% !!!
% Refrence t directions        [day]
t = (0:dt:(nt-1)*dt);
t_ref = repmat(t,nx,1);
t_up = t_ref/tmax;
t_step = t_up(1,2) - t_up(1,1);

% $$$
% Dimensionless x directions
x_less = x_ref/L;
x_step = 0.250;

% $$$
% Electric info
V = 25;                        % Voltage at z = 0 [V]
V_end = 0;                     % Voltage at z = L [V]
dVdx = (V-V_end)/L;            % Voltage gradient [V/m]

E_field = ones(nx,nt);         % Electric field [V]
M = linspace(V,0,nx);
for timestep = 1:nt
    E_field(:,timestep) = M;
end

% ###
% Setting inital value at t=0 and for all z (or x)
E_field(:,1) = (V-V_end);

% !!!
% Dimensionless Voltage     [Dimentionless]
phi_bar = E_field/(V-V_end);
dphidx = phi_bar./x_less;
dphidx(1,:) = 40; % Fixing damping value

% $$$
% Voltage gradient          [V/m]
E_field_dx = E_field/L;

% Global Physical info
T = 25 + 273;                  % Tempature [K]
R = 8.314;                     % Gas constant [J/mol.K]
F = 96485;                     % Faraady constant [C/mol]
D0 = 10^-9;                    % Reference diffusivity [m2/s]

% Soil info
n = 0.64;                      % Porosity   [Dimentionless]
tau = 1.25;                    % Tortuosity [Dimentionless]
dzdx = 1/tau;                  % divertion field

% Acetic acid info
sigma_surface = 0.0013;        % Surface conductivit [mhos/m]
K_a = 1.75*10^-6;              % dissociation constant [mol/m3]
K_H2O = 10^-8;                 % dissociation constant [(mol/m3)2]
K_b = 1.75*10^-6;              % dissociation constant [mol/m3]
mu_a = 0.001;                  % Solution viscosity [kg/m.s]
epsilon = 7*10^10;             % Electrical permittivity [F/m]
zeta = -0.0027;                % Zeta Potential [V]
zeta_0 = 2.6205e-23;           % Refrence Zeta Potential [V]

K_ads = 0.075;                 % Exprimental adsorbing constant
coeff = 1/(1+K_ads);

% Dimentionless     [Dimentionless]
Peclet = 47;
Beta = 967;
Z = 0.049;
k_0 = K_a;

% Valency           [Dimentionless]     $$$
z_HA = 0;
z_OH =-1;
z_Na = 1;
z_Cl =-1;
z_A = -1;
z_H = 1;
z_C = 0;

% $$$
% Species diffusivities (remapped)  [Dimentionless]
D_HA_less = 1.2;          % Acetic Acid
D_OH_less = 2.00;         % OH-
D_Na_less = 1.34;         % Na+
D_Cl_less = 2.05;         % Cl-
D_A_less = 1.2;           % Acid Agent
D_H_less = 9.35;          % H+
D_C_less = 2.09;          % Carbon

% $$$
% Species diffusivities (Normal)    [m2/s]      
D_HA = D_HA_less*D0;      % Acetic Acid
D_OH = D_OH_less*D0;      % OH-
D_Na = D_Na_less*D0;      % Na+
D_Cl = D_Cl_less*D0;      % Cl-
D_A = D_A_less*D0;        % Acid Agent
D_H = D_H_less*D0;        % H+
D_C = D_C_less*D0;        % Carbon

% $$$
% Species diffusivities (Normal)    [m2/day]      
D_HA_day = D_HA*24*3600;  % Acetic Acid
D_OH_day = D_OH*24*3600;  % OH-
D_Na_day = D_Na*24*3600;  % Na+
D_Cl_day = D_Cl*24*3600;  % Cl-
D_A_day = D_A*24*3600;    % Acid Agent
D_H_day = D_H*24*3600;    % H+
D_C_day = D_C*24*3600;    % Carbon

% $$$
% Mobility (Normal)       % [s·mol/kg]
v_HA = D_HA/(R*T);
v_OH = D_OH/(R*T);
v_Na = D_Na/(R*T);
v_Cl = D_Cl/(R*T);
v_A = D_A/(R*T);
v_H = D_H/(R*T);
v_C = D_C/(R*T);

% $$$
% Electroelectromigration velocity (Normal)     [m/s]
u_e_HA = -v_HA*z_HA*F*E_field_dx*(1/tau^2);
u_e_OH = -v_OH*z_OH*F*E_field_dx*(1/tau^2);
u_e_Na = -v_Na*z_Na*F*E_field_dx*(1/tau^2);
u_e_Cl = -v_Cl*z_Cl*F*E_field_dx*(1/tau^2);
u_e_A = -v_A*z_A*F*E_field_dx*(1/tau^2);
u_e_H = -v_H*z_H*F*E_field_dx*(1/tau^2);
u_e_C = -v_C*z_C*F*E_field_dx*(1/tau^2);

% $$$
% Dimentionless Calculated     [Dimentionless]
Peclet_calculated = (epsilon*zeta_0*V)/(mu_a*D0);
Beta_calculated = (F*V)/(R*T);
Z_calculated = (R*T*epsilon*zeta_0)/(D0*F*mu_a);

% $$$
% Refrence velocity            [m/s]
u_0 = (1/tau^2)*((epsilon*zeta_0)/mu_a)*((V-V_end)/L);

% $$$
% Convection velocity          [m/s]
u_eo = ((epsilon*zeta)/mu_a)*E_field_dx;
u_x = (epsilon/mu_a)*(zeta*E_field_dx);

% $$$
% Convection velocity (itself)
u_c = (1/tau^2)*(epsilon/mu_a)*(zeta*E_field_dx);
u_c1 = u_x/((tau^2));
% both are same and are here for assumptions

% $$$
% Convection velocity [Dimentionless]
u_c_up = -((zeta/zeta_0)*dphidx);
u_c_up1 = u_c/u_0;

% ###
% Testing new abberation
u_c = u_c*10^-17;
u_c_up1 = u_c_up1*10^-18;
u_c_up = u_c_up*10^-18;

% $$$
% Electroelectromigration velocity (1st)  [Dimentionless]
u_e_HA_less = (-1/Z_calculated)*D_HA_less*z_HA*dphidx;
u_e_OH_less = (-1/Z_calculated)*D_OH_less*z_OH*dphidx;
u_e_Na_less = (-1/Z_calculated)*D_Na_less*z_Na*dphidx;
u_e_Cl_less = (-1/Z_calculated)*D_Cl_less*z_Cl*dphidx;
u_e_A_less = (-1/Z_calculated)*D_A_less*z_A*dphidx;
u_e_H_less = (-1/Z_calculated)*D_H_less*z_H*dphidx;
u_e_C_less = (-1/Z_calculated)*D_C_less*z_C*dphidx;

% $$$
% Electroelectromigration velocity (2nd)  [Dimentionless]
u_e_HA_less1 = u_e_HA/u_0;
u_e_OH_less1 = u_e_OH/u_0;
u_e_Na_less1 = u_e_Na/u_0;
u_e_Cl_less1 = u_e_Cl/u_0;
u_e_A_less1 = u_e_A/u_0;
u_e_H_less1 = u_e_H/u_0;
u_e_C_less1 = u_e_C/u_0;
% These two need tests as they are different as result there are two total flux

% $$$
% Toatal velocity term (Normal)           [m/s]
u_t_HA = (u_e_HA + u_c);
u_t_OH = (u_e_OH + u_c);
u_t_Na = (u_e_Na + u_c);
u_t_Cl = (u_e_Cl + u_c);
u_t_H = (u_e_H + u_c);
u_t_A = (u_e_A + u_c);
u_t_C = (u_e_C + u_c);

% $$$
% Toatal velocity term (Normal)           [m/day]
u_t_HA_day = (u_e_HA + u_c)*24*3600;
u_t_OH_day = (u_e_OH + u_c)*24*3600;
u_t_Na_day = (u_e_Na + u_c)*24*3600;
u_t_Cl_day = (u_e_Cl + u_c)*24*3600;
u_t_H_day = (u_e_H + u_c)*24*3600;
u_t_A_day = (u_e_A + u_c)*24*3600;
u_t_C_day = (u_e_C + u_c)*24*3600;

% $$$
% Toatal velocity term (1st)
u_t_HA_up = (u_e_HA_less + u_c_up);
u_t_OH_up = (u_e_OH_less + u_c_up);
u_t_Na_up = (u_e_Na_less + u_c_up);
u_t_Cl_up = (u_e_Cl_less + u_c_up);
u_t_A_up = (u_e_A_less + u_c_up);
u_t_H_up = (u_e_H_less + u_c_up);
u_t_C_up = (u_e_C_less + u_c_up);

% $$$
% Toatal velocity term (2nd)
u_t_HA_up1 = (u_e_HA_less1 + u_c_up1);
u_t_OH_up1 = (u_e_OH_less1 + u_c_up1);
u_t_Na_up1 = (u_e_Na_less1 + u_c_up1);
u_t_Cl_up1 = (u_e_Cl_less1 + u_c_up1);
u_t_A_up1 = (u_e_A_less1 + u_c_up1);
u_t_H_up1 = (u_e_H_less1 + u_c_up1);
u_t_C_up1 = (u_e_C_less1 + u_c_up1);

% $$$
% Initial concentration        [mol/m3]
c_0 = 500; 
c_p = 200;
c_Na = c_p;
c_Cl = c_Na;

% $$$
% Initial Hydrocarbon concentration        [mg/kg]
c_C_TPH = 10000;

% $$$
% Hydrocarbon properties
API = 29.6;
MW = (6048/(API-5.9));                   % [g/mol]
rho = 1760;                              % [kg/(m3)]
bolian = 10^-3;                          % [g/mg]
c_C = ((c_C_TPH*rho*bolian)/MW);         % [mol/m3]


% !!!
% Sigma refrence               [S/m]
sigma_0 = ((F^2)*D0*c_0*(tau^2))/(R*T);

% !!!
% Alpha advection
alpha = u_0/(c_0*k_0);

% $$$
% --- Create arrays to save data for export
x_array = linspace(0,L,nx);
t_array = linspace(0,tmax,nt);

% $$$
% Erorr values for each finite step
h1 = (dt/(2*dx));
h2 = (dt/dx^2);

% --- Start of EKR (defualt)

% $$$
% Flux arrays
J_HA = zeros(nx,nt);
J_OH = zeros(nx,nt);
J_Na = zeros(nx,nt);
J_Cl = zeros(nx,nt);
J_A = zeros(nx,nt);
J_H = zeros(nx,nt);
J_C = zeros(nx,nt);

% $$$
% concentration arrays
G_HA = zeros(nx,nt);
G_OH = zeros(nx,nt);
G_Na = zeros(nx,nt);
G_Cl = zeros(nx,nt);
G_A = zeros(nx,nt);
G_H = zeros(nx,nt);
G_C = zeros(nx,nt);

% !!!
% concentration arrays [remapped]
G_HA_up = zeros(nx,nt);
G_OH_up = zeros(nx,nt);
G_Na_up = zeros(nx,nt);
G_Cl_up = zeros(nx,nt);
G_A_up = zeros(nx,nt);
G_H_up = zeros(nx,nt);
G_C_up = zeros(nx,nt);

% !!!
% adsorbed concentration arrays [remapped]
G_HA_ads = zeros(nx,nt);
G_OH_ads = zeros(nx,nt);
G_Na_ads = zeros(nx,nt);
G_Cl_ads = zeros(nx,nt);
G_A_ads = zeros(nx,nt);
G_H_ads = zeros(nx,nt);
G_C_ads = zeros(nx,nt);

% !!!
% Current array
i_z = zeros(nx, nt);

% !!!
% Sigma array
sigma_bar = zeros(nx,nt);
sum_HA = zeros(nx,nt);
sum_OH = zeros(nx,nt);
sum_Na = zeros(nx,nt);
sum_Cl = zeros(nx,nt);
sum_A = zeros(nx,nt);
sum_H = zeros(nx,nt);
sum_C = zeros(nx,nt);
sum_total = zeros(nx,nt);

% !!!
C_HA = zeros(nx,nt);
C_OH = zeros(nx,nt);
C_Na = zeros(nx,nt);
C_Cl = zeros(nx,nt);
C_A = zeros(nx,nt);
C_H = zeros(nx,nt);
C_C = zeros(nx,nt);
C_total = zeros(nx,nt);

% !!!
% SpeciesConstants
K_H2O_m = ones(nx,nt);
K_H2O_m = K_H2O_m*K_H2O;

K_a_m = ones(nx,nt);
K_a_m = K_a_m*K_a;

K_b_m = ones(nx,nt);
K_b_m = K_b_m*K_b;

% $$$
% Rate arrays
R_HA = zeros(nx,nt);
R_OH = zeros(nx,nt);
R_Na = zeros(nx,nt);
R_Cl = zeros(nx,nt);
R_A = zeros(nx,nt);
R_H = zeros(nx,nt);
R_C = zeros(nx,nt);

% !!!
% Adsorbed Rate arrays
R_HA_ads = zeros(nx,nt);
R_OH_ads = zeros(nx,nt);
R_Na_ads = zeros(nx,nt);
R_Cl_ads = zeros(nx,nt);
R_A_ads = zeros(nx,nt);
R_H_ads = zeros(nx,nt);
R_C_ads = zeros(nx,nt);

% $$$
% --- Set IC and BC
G_HA(:,1) = c_0;
G_OH(:,1) = c_0;
G_Na(:,1) = c_Na;
G_Cl(:,1) = c_Cl;
G_A(:,1) = c_0;
G_H(:,1) = c_0;
G_C(:,1) = c_C;

% !!!
G_HA_up(:,1) = G_HA(:,1)/c_0;
G_OH_up(:,1) = G_OH(:,1)/c_0;
G_Na_up(:,1) = G_Na(:,1)/c_Na;
G_Cl_up(:,1) = G_Cl(:,1)/c_Cl;
G_A_up(:,1) = G_A(:,1)/c_0;
G_H_up(:,1) = G_H(:,1)/c_0;
G_C_up(:,1) = G_C(:,1)/c_C;

% !!!
G_HA_ads(:,1) = K_ads*c_0;
G_OH_ads(:,1) = K_ads*c_0;
G_Na_ads(:,1) = K_ads*c_Na;
G_Cl_ads(:,1) = K_ads*c_Cl;
G_A_ads(:,1) = K_ads*c_0;
G_H_ads(:,1) = K_ads*c_0;
G_C_ads(:,1) = K_ads*c_C;

% !!!
% Fixing alpha advection
alpha_HA = (D_HA/(tau^2))*10^5;
alpha_OH = (D_OH/(tau^2))*10^5;
alpha_Na = (D_Na/(tau^2))*10^5;
alpha_Cl = (D_Cl/(tau^2))*10^5;
alpha_A = (D_A/(tau^2))*10^5;
alpha_H = (D_H/(tau^2))*10^5;
alpha_C = (D_C/(tau^2))*10^5;

% --- Set IC
G_HA(:,1) = c_0;
J_HA(1,:) = u_t_HA(1,:)*c_0;

G_OH(:,1) = c_0;
J_OH(1,:) = u_t_OH(1,:)*c_0;

G_Na(:,1) = c_0;
J_Na(1,:) = u_t_Na(1,:)*c_Na;

G_Cl(:,1) = c_0;
J_Cl(1,:) = u_t_Cl(1,:)*c_Cl;

G_A(:,1) = c_0;
J_A(1,:) = u_t_A(1,:)*c_0;

G_H(:,1) = c_0;
J_H(1,:) = u_t_H(1,:)*c_0;

G_C(:,1) = c_0;
J_C(1,:) = u_t_C(1,:)*c_C;

for m=1:nt-1
    % --- Set BC
    G_OH(1,m) = J_OH(1,m);   %--- Upper boundary
    G_HA(1,m) = J_HA(1,m);   %--- Upper boundary
    G_Na(1,m) = J_Na(1,m);   %--- Upper boundary
    G_Cl(1,m) = J_Cl(1,m);   %--- Upper boundary
    G_C(1,m) = J_C(1,m);     %--- Upper boundary
    G_H(1,m) = J_H(1,m);     %--- Upper boundary
    G_A(1,m) = J_A(1,m);     %--- Upper boundary

    for i=2:nx-1

        G_HA(i,m+1) = G_HA(i,m) + ((D_HA*h2)/tau^2)*(G_HA(i+1,m) -2*G_HA(i,m) + G_HA(i-1,m)) - h1*((u_t_HA(i+1,m) - u_t_HA(i,m))*(G_HA(i+1,m) - G_HA(i,m)));
        G_Na(i,m+1) = G_Na(i,m) + ((D_Na*h2)/tau^2)*(G_Na(i+1,m) -2*G_Na(i,m) + G_Na(i-1,m)) - h1*((u_t_Na(i+1,m) - u_t_Na(i,m))*(G_Na(i+1,m) - G_Na(i,m)));
        G_Cl(i,m+1) = G_Cl(i,m) + ((D_Cl*h2)/tau^2)*(G_Cl(i+1,m) -2*G_Cl(i,m) + G_Cl(i-1,m)) - h1*((u_t_Cl(i+1,m) - u_t_Cl(i,m))*(G_Cl(i+1,m) - G_Cl(i,m)));
        G_OH(i,m+1) = G_OH(i,m) + ((D_OH*h2)/tau^2)*(G_OH(i+1,m) -2*G_OH(i,m) + G_OH(i-1,m)) - h1*((u_t_OH(i+1,m) - u_t_OH(i,m))*(G_OH(i+1,m) - G_OH(i,m)));
        G_A(i,m+1) = G_A(i,m) + ((D_A*h2)/tau^2)*(G_A(i+1,m) -2*G_A(i,m) + G_A(i-1,m)) - h1*((u_t_A(i+1,m) - u_t_A(i,m))*(G_A(i+1,m) - G_A(i,m)));
       %G_H(i,m+1) = G_H(i,m) + ((D_H*h2)/tau^2)*(G_H(i+1,m) -2*G_H(i,m) + G_H(i-1,m)) - h1*((u_t_H(i+1,m) - u_t_H(i,m))*(G_H(i+1,m) - G_H(i,m)));
        G_C(i,m+1) = G_C(i,m) + ((D_C*h2)/tau^2)*(G_C(i+1,m) -2*G_C(i,m) + G_C(i-1,m)) - h1*((u_t_C(i+1,m) - u_t_C(i,m))*(G_C(i+1,m) - G_C(i,m)));
        
        G_H(i,m+1) = G_H(i,m) + ((D_H_day*h2)/tau^2)*(G_H(i+1,m) -2*G_H(i,m) + G_H(i-1,m)) - h1*((u_t_H(i,m))*(G_H(i+1,m) - G_H(i,m)));

        G_HA(end,m) = J_HA(i,m);    %--- Lower boundary
        G_OH(end,m) = J_OH(i,m);    %--- Lower boundary
        G_Na(end,m) = J_Na(i,m);    %--- Lower boundary
        G_Cl(end,m) = J_Cl(i,m);    %--- Lower boundary
        G_A(end,m) = J_A(i,m);      %--- Lower boundary
        G_H(end,m) = J_H(i,m);      %--- Lower boundary
        G_C(end,m) = J_C(i,m);      %--- Lower boundary


    
    end
end
% Convert concentration to mol/m3
G_HA_converted = G_HA_up*c_0;
G_A_converted = G_A_up*c_0;
G_Na_converted = G_Na_up*c_Na;
G_Cl_converted = G_Cl_up*c_Cl;
G_H_converted = G_H_up*c_0;
G_OH_converted = G_OH_up*c_0;
G_C_converted = G_C_up*c_C;

G_C_TPH_f = G_C_converted*(MW/(rho*bolian));

pH = log10(G_H_converted);

pH_scale = linspace(1,40,41);
xl = [0,5,10,15,20,25,30,35];
yl = [10000,7900,7100,6000,5700,5500,5400,5100];

figure(1);

hold on;


%plot(t_array,G_HA_converted(10,:),'-','DisplayName', 'HA');
%plot(t_array,G_A_converted(10,:),'-','DisplayName', 'A-');
%plot(t_array,G_Na_converted(10,:),'-','DisplayName', 'Na+');
%plot(t_array,G_Cl_converted(10,:),'-','DisplayName', 'Cl');
plot(t_array,G_H_converted(10,:),'-','DisplayName', 'H+ ');
%plot(t_array,G_OH_converted(10,:),'-','DisplayName', 'OH- ');
%plot(t_array,G_C_TPH_f(10,:),'-','DisplayName', 'Hydrocarbon');

%plot(t_array,R_H(10,:),'-','DisplayName', 'H+++ ');

legend();
