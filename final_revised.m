% Geomesh info
L = 0.4;                       % length of domain in x direction [m]
L_cm = L*100;                  % length of domain in x direction [cm]

% Nodes
tmax = 35;                     % end time [day]
nx = 41;                       % number of nodes in x direction
nt = 50401;                    % number of time steps
dx = L/(nx-1);                 % [m]
dt = tmax/(nt-1);              % [day]

% Refrence x directions        [m]
x = (10^-5:dx:(nx)*dx);
x = transpose(x);
x_ref = repmat(x,1,nt);

% Refrence t directions        [m]
t = (0:dt:(nt-1)*dt);
t_ref = repmat(t,nx,1);
t_up = t_ref/tmax;
t_step = t_up(1,2) - t_up(1,1);

% Electric info
phi = 25;                       % Voltage at 1st cap [V]
phi_end = 0;                    % Voltage at 2nd cap [V]
dEdx = phi/L;                   % Voltage gradient [V/m]

E_field = ones(nx,nt);          % Electric field [V]
M = linspace(phi,0,nx);
for timestep = 1:nt
    E_field(:,timestep) = M;
end

E_field_dx = E_field/L;         % Electric field in lenght [V/m]

% Global Physical info
T = 25 + 273;                   % Temperature [K]
R = 8.314;                      % Gas constant [J/mol.K]
F = 96485;                      % Faraady constant [C/mol]
D0 = 10^-9;                     % Reference diffusivity [m2/s]

% Soil info
n = 0.64;                       % Porosity
tau = 1.25;                     % Tortuosity
dzdx = 1/tau;                   % divertion field

% Acetic acid info
sigma_surface = 0.0013;         % Surface conductivit [mhos/m]
K_a = 1.75*10^-6;               % dissociation constant [mol/m3]
K_H2O = 10^-8;                  % dissociation constant [(mol/m3)2]
K_b = 1.75*10^-6;               % dissociation constant [mol/m3]
mu_a = 0.001;                   % Solution viscosity [kg/(m.day)]
epsilon = 7*10^10;              % Electrical permitivity [F/m]
zeta = -0.0027;                 % Zeta Potential [V]
zeta_0 = 2.6205e-23;            % Refrence Zeta Potential [V]

k_ads = 0.075;                  % Exprimental adsorbing constant
coeff = 1/(1+k_ads);            % Adsorbing coefficent
K = 0.02;                       % Exprimental Microbal constant
bolian = 24*3600;               % Time conversion cofactor

% Dimensionless parameters
Pe = 47;
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

z_total = z_HA + z_A + z_Na + z_Cl + z_H + z_OH + z_C;

% Species diffusivities (remapped)  [Dimentionless]
D_HA_Dless = 1.2;          % Acetic Acid
D_A_Dless = 1.2;           % Acid Agent
D_Na_Dless = 1.34;         % Na+
D_Cl_Dless = 2.05;         % Cl-
D_H_Dless = 9.35;          % H+
D_OH_Dless = 2.00;         % OH-
D_C_Dless = 2.00;          % Carbon

D_i_Dless = [D_HA_Dless, D_A_Dless, D_Na_Dless, D_Cl_Dless, D_H_Dless, D_OH_Dless, D_C_Dless];

% Species diffusivities (Normal)    [m2/s]
D_HA = D_HA_Dless*D0;      % Acetic Acid
D_A = D_A_Dless*D0;        % Acid Agent
D_Na = D_Na_Dless*D0;      % Na+
D_Cl = D_Cl_Dless*D0;      % Cl-
D_H = D_H_Dless*D0;        % H+
D_OH = D_OH_Dless*D0;      % OH-
D_C = D_C_Dless*D0;        % Carbon

D_i = [D_HA, D_A, D_Na, D_Cl, D_H, D_OH, D_C];

% Mobility (Normal)        % [s·mol/kg]
v_HA = D_HA/(R*T);
v_Na = D_Na/(R*T);
v_Cl = D_Cl/(R*T);
v_OH = D_OH/(R*T);
v_A = D_A/(R*T);
v_H = D_H/(R*T);
v_C = D_C/(R*T);

% Mobility (remapped)      % [s·mol/kg]
v_HA_Dless = D_HA_Dless/(R*T);
v_Na_Dless = D_Na_Dless/(R*T);
v_Cl_Dless = D_Cl_Dless/(R*T);
v_OH_Dless = D_OH_Dless/(R*T);
v_A_Dless = D_A_Dless/(R*T);
v_H_Dless = D_H_Dless/(R*T);
v_C_Dless = D_C_Dless/(R*T);

%-------------------------------
% Converting [s] into [day] is needed otherwise calculated D_standalone is wrong
% Species Diffuision abberation standalone      [Dimentionless]
D_standalone_OH = D_OH*bolian*(dt/dx^2);
D_standalone_HA = D_HA*bolian*(dt/dx^2);
D_standalone_Na = D_Na*bolian*(dt/dx^2);
D_standalone_Cl = D_Cl*bolian*(dt/dx^2);
D_standalone_A = D_A*bolian*(dt/dx^2);
D_standalone_H = D_H*bolian*(dt/dx^2);
D_standalone_C = D_C*bolian*(dt/dx^2);

% Species Diffusion advection standalone        [Dimentionless]
alpha_OH = D_standalone_OH/(n*tau^2);
alpha_HA = D_standalone_HA/(n*tau^2);
alpha_Na = D_standalone_Na/(n*tau^2);
alpha_Cl = D_standalone_Cl/(n*tau^2);
alpha_A = D_standalone_A/(n*tau^2);  
alpha_H = D_standalone_H/(n*tau^2);
alpha_C = D_standalone_C/(n*tau^2);

% Mobility (Normal)      % [s·mol/kg]
v_HA = D_HA/(R*T);
v_OH = D_OH/(R*T);
v_Na = D_Na/(R*T);
v_Cl = D_Cl/(R*T);
v_A = D_A/(R*T);
v_H = D_H/(R*T);
v_C = D_C/(R*T);

% Electroelectromigration velocity (Normal)     [m/s]
u_e_HA = -v_HA*z_HA*F*E_field*(1/tau^2);
u_e_OH = -v_OH*z_OH*F*E_field*(1/tau^2);
u_e_Na = -v_Na*z_Na*F*E_field*(1/tau^2);
u_e_Cl = -v_Cl*z_Cl*F*E_field*(1/tau^2);
u_e_A = -v_A*z_A*F*E_field*(1/tau^2);
u_e_H = -v_H*z_H*F*E_field*(1/tau^2);
u_e_C = -v_C*z_C*F*E_field*(1/tau^2);


% Refrence velocity                             [m/s]
u_0 = (1/tau^2)*((epsilon*zeta)/mu_a)*E_field;

% Convection velocity                           [m/s]
u_x = (epsilon/mu_a)*(zeta*E_field);

% Convection velocity (itself)
u_c = u_x/((tau^2)*10^19);

% Toatal velocity term (Normal)
u_t_HA = (u_e_HA + u_c);
u_t_OH = (u_e_OH + u_c);
u_t_Na = (u_e_Na + u_c);
u_t_Cl = (u_e_Cl + u_c);
u_t_A = (u_e_A + u_c);
u_t_H = (u_e_H + u_c);
u_t_C = (u_e_C + u_c);

% Velocity advection without coefficent
beta_OH = u_t_OH*(dt/2*dx);
beta_HA = u_t_HA*(dt/2*dx);
beta_Na = u_t_Na*(dt/2*dx);
beta_Cl = u_t_Cl*(dt/2*dx);
beta_C = u_t_C*(dt/2*dx);
beta_H = u_t_H*(dt/2*dx);
beta_A = u_t_A*(dt/2*dx);


% Initial concentration        [mol/m3]
c_0 = 500; 
c_p = 200;
c_Na = c_p;
c_Cl = c_Na;

% Initial Hydrocarbon concentration        [mg/kg]
c_C_TPH = 10000;

% Hydrocarbon properties
API = 29.6;
MW = (6048/(API-5.9));                   % [g/mol]
rho = 1760;                              % [kg/(m3)]
bolian = 10^-3;                          % [g/mg]
c_C = ((c_C_TPH*rho*bolian)/MW);         % [mol/m3]

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
G_C(:,1) = c_C;

G_HA_up(:,1) = G_HA(:,1)/c_0;
G_A_up(:,1) = G_A(:,1)/c_0;
G_Na_up(:,1) = G_Na(:,1)/c_Na;
G_Cl_up(:,1) = G_Cl(:,1)/c_Cl;
G_H_up(:,1) = G_H(:,1)/c_0;
G_OH_up(:,1) = G_OH(:,1)/c_0;
G_C_up(:,1) = G_C(:,1)/c_C;

G_HA_ads(:,1) = K_ads*c_0;
G_A_ads(:,1) = K_ads*c_0;
G_Na_ads(:,1) = K_ads*c_Na;
G_Cl_ads(:,1) = K_ads*c_Cl;
G_H_ads(:,1) = K_ads*c_0;
G_OH_ads(:,1) = K_ads*c_0;
G_C_ads(:,1) = K_ads*c_C;

% Fixing alpha advection
alpha_HA = (D_HA/(tau^2))*10^5;
alpha_A = (D_A/(tau^2))*10^5;
alpha_Na = (D_Na/(tau^2))*10^5;
alpha_Cl = (D_Cl/(tau^2))*10^5;
alpha_H = (D_H/(tau^2))*10^5;
alpha_OH = (D_OH/(tau^2))*10^5;
alpha_C = (D_C/(tau^2))*10^5;

% --- Set IC and BC

G_C(:,1)= 10000;
J_C(1,:)= (u_c_ekr(1,:) + u_e_C(1,:))*10000;

G_H(:,1)= 10000;
J_H(1,:)= (u_c_ekr(1,:) + u_e_H(1,:))*10000;

G_OH(:,1)= 10000;
J_OH(1,:)= (u_c_ekr(1,:) + u_e_OH(1,:))*10000;

G_HA(:,1)= 10000;
J_HA(1,:)= (u_c_ekr(1,:) + u_e_HA(1,:))*10000;

G_BOH(:,1)= 10000;
J_BOH(1,:)= (u_c_ekr(1,:) + u_e_BOH(1,:))*10000;

G_A(:,1)= 10000;
J_A(1,:)= (u_c_ekr(1,:) + u_e_A(1,:))*10000;

G_B(:,1)= 10000;
J_B(1,:)= (u_c_ekr(1,:) + u_e_B(1,:))*10000;

R_C(:,1) = R_i*(dt)/n;
R_OH(:,1) = R_i*(dt)/n;
R_H(:,1) = R_i*(dt)/n;
R_HA(:,1) = R_i*(dt)/n;
R_BOH(:,1) = R_i*(dt)/n;
R_B(:,1) = R_i*(dt)/n;
R_A(:,1) = R_i*(dt)/n;

s_H(:,1) = (z_H^2)*v_H*G_H(:,1);
s_OH(:,1) = (z_OH^2)*v_OH*G_OH(:,1);
s_C(:,1) = (z_C^2)*v_C*G_C(:,1);
Sigma(:,1) = (F^2)*(s_H(:,1) + s_OH(:,1) + s_C(:,1));

for m= 1:nt-1

    G_C(1,m) =J_C(1,m); %--- Upper boundary
    G_H(1,m) =J_H(1,m); %--- Upper boundary
    G_OH(1,m) =J_OH(1,m); %--- Upper boundary
    G_HA(1,m) =J_HA(1,m); %--- Upper boundary
    G_BOH(1,m) =J_BOH(1,m); %--- Upper boundary
    G_A(1,m) =J_A(1,m); %--- Upper boundary
    G_B(1,m) =J_B(1,m); %--- Upper boundary

    for i= 2:nx-1
        
        G_C(i,m+1) = G_C(i,m) + alpha_C*(G_C(i+1,m) -2*G_C(i,m) + G_C(i-1,m)) + beta_C(i,m)*(G_C(i+1,m) - G_C(i-1,m)) + R_C(i,m)/R_D;
        G_H(i,m+1) = G_H(i,m) + alpha_H*(G_H(i+1,m) -2*G_H(i,m) + G_H(i-1,m)) + beta_H(i,m)*(G_H(i+1,m) - G_H(i-1,m)) + R_H(i,m)/R_D;
        G_OH(i,m+1) = G_OH(i,m) + alpha_OH*(G_OH(i+1,m) -2*G_OH(i,m) + G_OH(i-1,m)) + beta_OH(i,m)*(G_OH(i+1,m) - G_OH(i-1,m)) + R_OH(i,m)/R_D;
        G_HA(i,m+1) = G_HA(i,m) + alpha_HA*(G_HA(i+1,m) -2*G_HA(i,m) + G_HA(i-1,m)) + beta_HA(i,m)*(G_HA(i+1,m) - G_HA(i-1,m)) + R_HA(i,m)/R_D;
        G_BOH(i,m+1) = G_BOH(i,m) + alpha_BOH*(G_BOH(i+1,m) -2*G_BOH(i,m) + G_BOH(i-1,m)) + beta_BOH(i,m)*(G_BOH(i+1,m) - G_BOH(i-1,m)) + R_BOH(i,m)/R_D;
        G_A(i,m+1) = G_A(i,m) + alpha_A*(G_A(i+1,m) -2*G_A(i,m) + G_A(i-1,m)) + beta_A(i,m)*(G_A(i+1,m) - G_A(i-1,m)) + R_A(i,m)/R_D;
        G_B(i,m+1) = G_B(i,m) + alpha_B*(G_B(i+1,m) -2*G_B(i,m) + G_B(i-1,m)) + beta_B(i,m)*(G_B(i+1,m) - G_B(i-1,m)) + R_B(i,m)/R_D;

        J_C(i,m) = (u_c_ekr(i,m) + u_e_C(i,m))*G_C(i,m) - alpha_C*(G_C(i+1,m) - G_C(i,m));
        J_H(i,m) = (u_c_ekr(i,m) + u_e_H(i,m))*G_H(i,m) - alpha_H*(G_H(i+1,m) - G_H(i,m));
        J_OH(i,m) = (u_c_ekr(i,m) + u_e_OH(i,m))*G_OH(i,m) - alpha_OH*(G_OH(i,m) - G_OH(i-1,m));
        J_HA(i,m) = (u_c_ekr(i,m) + u_e_HA(i,m))*G_HA(i,m) - alpha_HA*(G_HA(i,m) - G_HA(i-1,m));
        J_BOH(i,m) = (u_c_ekr(i,m) + u_e_BOH(i,m))*G_BOH(i,m) - alpha_BOH*(G_BOH(i,m) - G_BOH(i-1,m));
        J_A(i,m) = (u_c_ekr(i,m) + u_e_A(i,m))*G_A(i,m) - alpha_A*(G_A(i+1,m) - G_A(i,m));
        J_B(i,m) = (u_c_ekr(i,m) + u_e_B(i,m))*G_B(i,m) - alpha_B*(G_B(i+1,m) - G_B(i,m));
        
        G_C(end,m) = J_C(i,m); %--- Lower boundary
        G_H(end,m) = J_H(i,m); %--- Lower boundary
        G_OH(end,m) = J_OH(i,m); %--- Lower boundary
        G_HA(end,m) = J_HA(i,m); %--- Lower boundary
        G_BOH(end,m) = J_BOH(i,m); %--- Lower boundary
        G_A(end,m) = J_A(i,m); %--- Lower boundary
        G_B(end,m) = J_B(i,m); %--- Lower boundary

        s_H(i,m) = (z_H^2)*v_H*G_H(i,m);
        s_OH(i,m) = (z_OH^2)*v_OH*G_OH(i,m);
        s_C(i,m) = (z_C^2)*v_H*G_C(i,m);
        Sigma(i,m) = (F^2)*(s_C(i,m) + s_H(i,m) + s_OH(i,m)) + Sigma_ref(i,m);
        
        i_z(i,m) = (-1*Sigma(i,m)*E_field(i,m) - F*((z_C*D_C*(G_C(i+1,m) - G_C(i-1,m))) + (z_H*D_H*(G_H(i+1,m) - G_H(i-1,m))) + (z_OH*D_OH*(G_OH(i+1,m) - G_OH(i-1,m)))))/(tau^2);
        if i == 2
            R_prime_H = i_z(i,m)/F;
            R_H(i,m) = -1*R_prime_H;
        end
        if i == nx-1
            R_prime_OH = i_z(i,m)/F;
            R_OH(i,m) = -1*R_prime_OH;
        end
        K_H2O(i,m) = G_H(i,m)*G_OH(i,m);
        K_a(i,m) = (G_H(i,m)*G_A(i,m))/G_HA(i,m);
        K_b(i,m) = (G_B(i,m)*G_OH(i,m))/G_BOH(i,m);
        R_H(i,m) = (K_H2O(i,m)*G_H(i,m)) + (K_a(i,m)*G_HA(i,m));
        R_OH(i,m) = (K_H2O(i,m)*G_OH(i,m)) + (K_b(i,m)*G_BOH(i,m));
        R_B(i,m) = (K_b(i,m)*G_BOH(i,m));
        R_A(i,m) = (K_a(i,m)*G_HA(i,m));
        R_C(i,m) = R_i*(dt)/n;
        
    end
end


G_C_B = zeros(nx,nt);
G_H_B = zeros(nx,nt);
G_OH_B = zeros(nx,nt);
G_HA_B = zeros(nx,nt);
G_BOH_B = zeros(nx,nt);
G_A_B = zeros(nx,nt);
G_B_B = zeros(nx,nt);

J_C_B = zeros(nx,nt);
J_H_B = zeros(nx,nt);
J_OH_B = zeros(nx,nt);
J_HA_B = zeros(nx,nt);
J_BOH_B = zeros(nx,nt);
J_A_B = zeros(nx,nt);
J_B_B = zeros(nx,nt);

K_H2O_B = zeros(nx,nt);
K_a_B = zeros(nx,nt);
K_b_B = zeros(nx,nt);

R_C_B = zeros(nx,nt);
R_H_B = zeros(nx,nt);
R_OH_B = zeros(nx,nt);
R_HA_B = zeros(nx,nt);
R_BOH_B = zeros(nx,nt);
R_A_B = zeros(nx,nt);
R_B_B = zeros(nx,nt);

Sigma_B = zeros(nx,nt);
Sigma_ref_B = ones(nx,nt);
sigma_ref_B = Sigma_ref_B*sigma_surface;

i_z_B = zeros(nx, nt);

s_H_B = zeros(nx,nt);
s_C_B = zeros(nx,nt);
s_OH_B = zeros(nx,nt);

% --- Set IC and BC

G_C_B(:,1)= 10000;
J_C_B(1,:)= (u_c(1,:) + u_e_C(1,:))*10000;

G_H_B(:,1)= 10000;
J_H_B(1,:)= (u_c(1,:) + u_e_H(1,:))*10000;

G_OH_B(:,1)= 10000;
J_OH_B(1,:)= (u_c(1,:) + u_e_OH(1,:))*10000;

G_HA_B(:,1)= 10000;
J_HA_B(1,:)= (u_c(1,:) + u_e_HA(1,:))*10000;

G_BOH_B(:,1)= 10000;
J_BOH_B(1,:)= (u_c(1,:) + u_e_BOH(1,:))*10000;

G_A_B(:,1)= 10000;
J_A_B(1,:)= (u_c(1,:) + u_e_A(1,:))*10000;

G_B_B(:,1)= 10000;
G_B_B(1,:)= (u_c(1,:) + u_e_B(1,:))*10000;

R_C_B(:,1) = R_i*coeff*(dt)/n;
R_OH_B(:,1) = R_i*coeff*(dt)/n;
R_H_B(:,1) = R_i*coeff*(dt)/n;
R_HA_B(:,1) = R_i*coeff*(dt)/n;
R_BOH_B(:,1) = R_i*coeff*(dt)/n;
R_B_B(:,1) = R_i*coeff*(dt)/n;
R_A_B(:,1) = R_i*coeff*(dt)/n;

s_H_B(:,1) = (z_H^2)*v_H*G_H_B(:,1);
s_OH_B(:,1) = (z_OH^2)*v_OH*G_OH_B(:,1);
s_C_B(:,1) = (z_C^2)*v_C*G_C_B(:,1);
Sigma_B(:,1) = (F^2)*(s_H_B(:,1) + s_OH_B(:,1) + s_C_B(:,1));

for xx = 1:nx
    for tt = 1:nt
        M_g = exp(-K*(tt/1440)*(Z*(Beta/Pe)));
        sub(xx,tt) = M_g;
        fixup(xx,tt) = 1/(-0.01*(tt/1440) + 0.61);
        %fixup(xx,tt) = 1/(-0.01*(tt/1440) + 0.61);
    end
end
growth = sub.*fixup;

for m= 1:nt-1

    G_C_B(1,m) =J_C_B(1,m); %--- Upper boundary
    G_H_B(1,m) =J_H_B(1,m); %--- Upper boundary
    G_OH_B(1,m) =J_OH_B(1,m); %--- Upper boundary
    G_HA_B(1,m) =J_HA_B(1,m); %--- Upper boundary
    G_BOH_B(1,m) =J_BOH_B(1,m); %--- Upper boundary
    G_A_B(1,m) =J_A_B(1,m); %--- Upper boundary
    G_B_B(1,m) = J_B_B(1,m); %--- Upper boundary

    for i= 2:nx-1
        
        
        G_C_B(i,m+1) = G_C_B(i,m) + growth(i,m)*(alpha_prime_C*(G_C_B(i+1,m) -2*G_C_B(i,m) + G_C_B(i-1,m)) + beta_prime_C(i,m)*(G_C_B(i+1,m) - G_C_B(i-1,m)) + R_C_B(i,m)/R_D);
        G_H_B(i,m+1) = G_H_B(i,m) + growth(i,m)*(alpha_prime_H*(G_H_B(i+1,m) -2*G_H_B(i,m) + G_H_B(i-1,m)) + beta_prime_H(i,m)*(G_H_B(i+1,m) - G_H_B(i-1,m)) + R_H_B(i,m)/R_D);
        G_OH_B(i,m+1) = G_OH_B(i,m) + growth(i,m)*(alpha_prime_OH*(G_OH_B(i+1,m) -2*G_OH_B(i,m) + G_OH_B(i-1,m)) + beta_prime_OH(i,m)*(G_OH_B(i+1,m) - G_OH_B(i-1,m)) + R_OH_B(i,m)/R_D);
        G_HA_B(i,m+1) = G_HA_B(i,m) + growth(i,m)*(alpha_prime_HA*(G_HA_B(i+1,m) -2*G_HA_B(i,m) + G_HA_B(i-1,m)) + beta_prime_HA(i,m)*(G_HA_B(i+1,m) - G_HA_B(i-1,m)) + R_HA_B(i,m)/R_D);
        G_BOH_B(i,m+1) = G_BOH_B(i,m) + growth(i,m)*(alpha_prime_BOH*(G_BOH_B(i+1,m) -2*G_BOH_B(i,m) + G_BOH_B(i-1,m)) + beta_prime_BOH(i,m)*(G_BOH_B(i+1,m) - G_BOH_B(i-1,m)) + R_BOH_B(i,m)/R_D);
        G_A_B(i,m+1) = G_A_B(i,m) + growth(i,m)*(alpha_prime_A*(G_A_B(i+1,m) -2*G_A_B(i,m) + G_A_B(i-1,m)) + beta_prime_A(i,m)*(G_A_B(i+1,m) - G_A_B(i-1,m)) + R_A_B(i,m)/R_D);
        G_B_B(i,m+1) = G_B_B(i,m) + growth(i,m)*(alpha_prime_B*(G_B_B(i+1,m) -2*G_B_B(i,m) + G_B_B(i-1,m)) + beta_prime_B(i,m)*(G_B_B(i+1,m) - G_B_B(i-1,m)) + R_B_B(i,m)/R_D);
        
        J_C_B(i,m) = (u_c(i,m) + u_e_C(i,m))*G_C_B(i,m) - alpha_prime_C*(G_C_B(i+1,m) - G_C_B(i,m));
        J_H_B(i,m) = (u_c(i,m) + u_e_H(i,m))*G_H_B(i,m) - alpha_prime_H*(G_H_B(i+1,m) - G_H_B(i,m));
        J_OH_B(i,m) = (u_c(i,m) + u_e_OH(i,m))*G_OH_B(i,m) - alpha_prime_OH*(G_OH_B(i,m) - G_OH_B(i-1,m));
        J_HA_B(i,m) = (u_c(i,m) + u_e_HA(i,m))*G_HA_B(i,m) - alpha_prime_HA*(G_HA_B(i,m) - G_HA_B(i-1,m));
        J_BOH_B(i,m) = (u_c(i,m) + u_e_BOH(i,m))*G_BOH_B(i,m) - alpha_prime_BOH*(G_BOH_B(i,m) - G_BOH_B(i-1,m));
        J_A_B(i,m) = (u_c(i,m) + u_e_A(i,m))*G_A_B(i,m) - alpha_prime_A*(G_A_B(i,m) - G_A_B(i-1,m));
        J_B_B(i,m) = (u_c(i,m) + u_e_B(i,m))*G_B_B(i,m) - alpha_prime_B*(G_B_B(i,m) - G_B_B(i-1,m));

        G_C_B(end,m) = J_C_B(i,m); %--- Lower boundary
        G_H_B(end,m) = J_H_B(i,m); %--- Lower boundary
        G_OH_B(end,m) = J_OH_B(i,m); %--- Lower boundary
        G_HA_B(end,m) = J_HA_B(i,m); %--- Lower boundary
        G_BOH_B(end,m) = J_BOH_B(i,m); %--- Lower boundary
        G_A_B(end,m) = J_A_B(i,m); %--- Lower boundary
        G_B_B(end,m) = J_B_B(i,m); %--- Lower boundary

        s_H_B(i,m) = (z_H^2)*v_H*G_H_B(i,m);
        s_OH_B(i,m) = (z_OH^2)*v_OH*G_OH_B(i,m);
        s_C_B(i,m) = (z_C^2)*v_H*G_C_B(i,m);
        Sigma_B(i,m) = (F^2)*(s_C_B(i,m) + s_H_B(i,m) + s_OH_B(i,m)) + Sigma_ref(i,m);
        
        i_z_B(i,m) = (-1*Sigma_B(i,m)*E_field(i,m) - F*((z_C*D_C*(G_C_B(i+1,m) - G_C(i-1,m))) + (z_H*D_H*(G_H_B(i+1,m) - G_H_B(i-1,m))) + (z_OH*D_OH*(G_OH_B(i+1,m) - G_OH_B(i-1,m)))))/(tau^2);
        if i == 2
            R_prime_H_B = i_z_B(i,m)/F;
            R_H_B(i,m) = -1*R_prime_H_B;
        end
        if i == nx-1
            R_prime_OH_B = i_z_B(i,m)/F;
            R_OH_B(i,m) = -1*R_prime_OH_B;
        end
        K_H2O_B(i,m) = G_H_B(i,m)*G_OH_B(i,m);
        K_a_B(i,m) = (G_H_B(i,m)*G_A_B(i,m))/G_HA_B(i,m);
        K_b_B(i,m) = (G_B_B(i,m)*G_OH_B(i,m))/G_BOH_B(i,m);
        R_H_B(i,m) = (K_H2O(i,m)*G_H_B(i,m)) + (K_a(i,m)*G_HA_B(i,m));
        R_OH_B(i,m) = (K_H2O(i,m)*G_OH_B(i,m)) + (K_b(i,m)*G_BOH_B(i,m));
        R_B_B(i,m) = (K_b(i,m)*G_BOH_B(i,m));
        R_A_B(i,m) = (K_a(i,m)*G_HA_B(i,m));
        R_C_B(i,m) = R_i*coeff*(dt)/n;
        
    end
end


pH = log10(G_H);
pH_B = log10(G_H_B);
x_scale = linspace(1,40,41);
xl = [0,5,10,15,20,25,30,35];
yl = [10000,7900,7100,6000,5700,5500,5400,5100];

figure(1);
hold on;

plot(t,G_C(10,:),'-','DisplayName', 'Hydrocarbon (EKR)');
plot(t,G_C_B(10,:),'--','DisplayName', 'Hydrocarbon (BKR)');

scatter(xl,yl, 'DisplayName', 'Expriment Data');

xlabel('Time');
ylabel('Conc(mg/kg)');

legend();
hold off;

% Plot for figure 2
gif_filename = 'pH_change_animation.gif';

figure(2);
hold on;

fig1 = plot(x_scale, pH(:,1), 'DisplayName', 'pH change in EKR');
fig2 = plot(x_scale, pH_B(:,1), 'DisplayName', 'pH change in BKR');

xlabel('Distance (cm)');
ylabel('pH Level');

for h = 1:50401
    set(fig1, 'XData', x_scale, 'YData', pH(:,h)');
    set(fig2, 'XData', x_scale, 'YData', pH_B(:,h)');

    h2 = h / 1440;
    title(['pH Change - Day: ', sprintf('%.2f', h2)]);
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if h == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
    end
    
    pause(0.00001);
end

legend;
hold off;



% Plot for figure 3
figure(3)
plot(x_scale,E_field(:,50400),'--','DisplayName', 'dEdx')
xlabel('Length (cm)');
ylabel('E field');
title('Electric feild gradient')
legend();