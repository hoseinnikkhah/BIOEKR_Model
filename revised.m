% Geomesh info
L = 0.4;                        % length of domain in x direction [m]       
tmax = 35;                      % end time [day]
nx = 41;                        % number of nodes in x direction
nt = 50401;                     % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);


% Physical info
T = 25 + 273;                   % Temperature [K]
R = 8.314;                      % Gas constant [J/mol.K]
n = 0.64;                       % Porosity
F = 96485;                      % Faraady constant [C/mol]
phi = 25;                       % Electric field [V]
phi_end = 0;
tau = 1.25;                     % Tortuosity
dzdx = 1/tau;                   % divertion field
dEdx = phi/L;                   % Electric field in lenght [V/m]
epsilon = 7*10^10;              % Electrical permitivity [F/m]
mu_oil = 510*24*3600;           % Oil viscosity [kg/(m.day)]
mu_solution = 0.001*24*3600;    % Solution viscosity [kg/(m.day)]
zeta = -0.0027;                 % Zeta potential [V]
K = 0.02;                       % Exprimental Microbal constant
K_A = 1.75*10^-5;               % Dissociation constant [mol/m3]
k_i = 0.075;                    % Exprimental constant
coeff = 1/(1+k_i);              % Adsorbing coefficent
R_i = (0.693/53.2);             % Initial Reaction flow rate

% Dimensionless parameters
Pe = 47;
Z = 0.049;
Beta = 967;

E_field = ones(nx,nt);
M = linspace(dEdx,0,nx);
for timestep = 1:nt
    E_field(:,timestep) = M;
end

E_field_dx = E_field/L;

% Species Valency
z_H = 1;
z_OH = -1;
z_C = 0;
z_HA = 0;
z_Na = 0;
z_A = -1;
z_Cl = 1;

% Current density
i = F*(z_H+z_OH+z_C+z_HA+z_Na+z_A+z_Cl);
% there is a flux term at the end as well but it is not calculated since sum of valencies are zero in this phenomena


% Counductivity
sigma_surface = 2.74*10^7;        % Conductivity [S/m]

% Reference diffusivity             [m2/s]
D0 = 10^-9;

% Species diffusivities (remapped)  [Dimentionless]
D_HA_less = 1.2;          % Acetic Acid
D_OH_less = 2.00;         % OH-
D_Na_less = 1.34;         % Na+
D_Cl_less = 2.05;         % Cl-
D_A_less = 1.2;           % Acid Agent
D_H_less = 9.35;          % H+
D_C_less = 2.00;          % Carbon


% Species diffusivities (Normal)    [m2/s]
D_HA = D_HA_less*D0;      % Acetic Acid
D_OH = D_OH_less*D0;      % OH-
D_Na = D_Na_less*D0;      % Na+
D_Cl = D_Cl_less*D0;      % Cl-
D_A = D_A_less*D0;        % Acid Agent
D_H = D_H_less*D0;        % H+
D_C = D_C_less*D0;        % Carbon

% Species diffusivities (Normal)    [m^2/day]
D_HA = D_HA*24*3600;      % Acetic Acid
D_OH = D_OH*24*3600;      % OH-
D_Na = D_Na*24*3600;      % Na+
D_Cl = D_Cl*24*3600;      % Cl-
D_A = D_A*24*3600;        % Acid Agent
D_H = D_H*24*3600;        % H+
D_C = D_C*24*3600;        % Carbon

% Species Diffuision abberation without coeff   [Dimentionless]
D_star_HA = D_HA*(dt/dx^2);
D_star_OH = D_OH*(dt/dx^2);
D_star_Na = D_Na*(dt/dx^2);
D_star_Cl = D_Cl*(dt/dx^2);
D_star_H = D_H*(dt/dx^2);
D_star_C = D_C*(dt/dx^2);
D_star_A = D_A*(dt/dx^2);
    
% Species Diffuision abberation with coeff      [Dimentionless]
D_prime_HA = D_HA*coeff*(dt/dx^2);
D_prime_OH = D_OH*coeff*(dt/dx^2);
D_prime_Na = D_Na*coeff*(dt/dx^2);
D_prime_Cl = D_Cl*coeff*(dt/dx^2);
D_prime_H = D_H*coeff*(dt/dx^2);
D_prime_C = D_C*coeff*(dt/dx^2);
D_prime_A = D_A*coeff*(dt/dx^2);

% Species Diffusion advection without coeff     [Dimentionless]
alpha_OH = D_star_OH/(n*tau^2);
alpha_HA = D_star_HA/(n*tau^2);
alpha_Na = D_star_Na/(n*tau^2);
alpha_Cl = D_star_Cl/(n*tau^2);
alpha_H = D_star_H/(n*tau^2);
alpha_C = D_star_C/(n*tau^2);
alpha_A = D_star_A/(n*tau^2);

% Species Diffusion advection with coeff       [Dimentionless]
alpha_prime_H = D_prime_H/(n*tau^2);
alpha_prime_C = D_prime_C/(n*tau^2);
alpha_prime_OH = D_prime_OH/(n*tau^2);
alpha_prime_HA = D_prime_HA/(n*tau^2);
alpha_prime_Na = D_prime_Na/(n*tau^2);
alpha_prime_A = D_prime_A/(n*tau^2);
alpha_prime_Cl = D_prime_Cl/(n*tau^2);

% Species Mobility                      [s·mol/kg]
v_OH = (D_OH/(R*T));
v_HA = (D_HA/(R*T));
v_Na = (D_Na/(R*T));
v_Cl = (D_Cl/(R*T));
v_C = (D_C/(R*T));
v_H = (D_H/(R*T));
v_A = (D_A/(R*T));

% Species electromigration velocity     [m/s] 
u_e_OH = (v_OH*z_OH*F*E_field_dx)/(tau^2);
u_e_HA = (v_HA*z_HA*F*E_field_dx)/(tau^2);
u_e_Na = (v_Na*z_Na*F*E_field_dx)/(tau^2);
u_e_Cl = (v_Cl*z_Cl*F*E_field_dx)/(tau^2);
u_e_H = (v_H*z_H*F*E_field_dx)/(tau^2); 
u_e_C = (v_C*z_C*F*E_field_dx)/(tau^2);
u_e_A = (v_A*z_A*F*E_field_dx)/(tau^2);
% note: E_field_dx is needed as it is dEdx in main formula

% domain velocity
u_x = (epsilon/mu_solution)*(zeta*E_field);     % Volumetric Velocity [m/s]

u_c = ones(nx,nt);
u_C = ((1/tau^2)*u_x)*Z/(24*3600*Pe*Beta);      % Convection Velocity [m/s]
u_C_ekr = ((1/tau^2)*u_x)*Z/(24*3600*Pe*Beta);  % Convection Velocity [m/s]

u_c = u_c.*u_C;
u_c_ekr = u_c.*u_C_ekr;

u_eo = ((epsilon*zeta)/mu_solution)*E_field;    % Electoosmotic Velocity [m3/s]
u_s = n*u_c;

e_r = 7.5;                                      % relative permittivity of clay [F/m]
e_0 = 8.854*10^-12;                             % permittivity of free space [F/m]

k_eo = (epsilon*zeta)/(mu*epsilon_oil);         % Electroosmotic mobility [m^2/V.s] Based on epsilon of both crude oil and clay, might be true form
k_eo1 = (e_0*zeta*recip)/(mu*(1-n)*(n^3));      % Electroosmotic mobility [m^2/V.s] Based on debye
k_eo2 = ((-7*10^-10)*zeta*n)/(mu*(tau^2));      % Electroosmotic mobility [m^2/V.s] Based on a journal (exact one)

% Species total velocity
u_t_OH = (u_e_OH + u_c)/n;
u_t_HA = (u_e_HA + u_c)/n;
u_t_Na = (u_e_Na + u_c)/n;
u_t_Cl = (u_e_Cl + u_c)/n;
u_t_C = (u_e_C + u_c)/n;
u_t_H = (u_e_H + u_c)/n;
u_t_A = (u_e_A + u_c)/n;


% Species total velocity
u_t_OH_ekr = (u_e_OH + u_c_ekr)/n;
u_t_HA_ekr = (u_e_HA + u_c_ekr)/n;
u_t_Na_ekr = (u_e_Na + u_c_ekr)/n;
u_t_Cl_ekr = (u_e_Cl + u_c_ekr)/n;
u_t_C_ekr = (u_e_C + u_c_ekr)/n;
u_t_H_ekr = (u_e_H + u_c_ekr)/n;
u_t_A_ekr = (u_e_A + u_c_ekr)/n;

% Velocity advection without coefficent     [Dimentionless]
beta_OH = u_t_OH_ekr*(dt/2*dx);
beta_HA = u_t_HA_ekr*(dt/2*dx);
beta_Na = u_t_Na_ekr*(dt/2*dx);
beta_Cl = u_t_Cl_ekr*(dt/2*dx);
beta_C = u_t_C_ekr*(dt/2*dx);
beta_H = u_t_H_ekr*(dt/2*dx);
beta_A = u_t_A_ekr*(dt/2*dx);

% Velocity advection with coefficent        [Dimentionless]
beta_prime_OH = coeff*u_t_OH*(dt/2*dx);
beta_prime_HA = coeff*u_t_HA*(dt/2*dx);
beta_prime_Na = coeff*u_t_Na*(dt/2*dx);
beta_prime_Cl = coeff*u_t_Cl*(dt/2*dx);
beta_prime_C = coeff*u_t_C*(dt/2*dx); 
beta_prime_H = coeff*u_t_H*(dt/2*dx);
beta_prime_A = coeff*u_t_A*(dt/2*dx);

R_D = coeff*(dt)/n;              % Reaction rate Dimensionless factor

% --- Create arrays to save data for export
x_array = linspace(0,L,nx);
t_array = linspace(0,tmax,nt);

J_OH = zeros(nx,nt);
J_HA = zeros(nx,nt);
J_Na = zeros(nx,nt);
J_Cl = zeros(nx,nt);
J_C = zeros(nx,nt);
J_H = zeros(nx,nt);
J_A = zeros(nx,nt);

G_OH = zeros(nx,nt);
G_HA = zeros(nx,nt);
G_Na = zeros(nx,nt);
G_Cl = zeros(nx,nt);
G_C = zeros(nx,nt);
G_H = zeros(nx,nt);
G_A = zeros(nx,nt);

Sigma = zeros(nx,nt);
Sigma_ref = ones(nx,nt);
sigma_ref = Sigma_ref*sigma_surface;

i_z = zeros(nx, nt);

s_H = zeros(nx,nt);
s_C = zeros(nx,nt);
s_OH = zeros(nx,nt);

K_H2O = zeros(nx,nt);
K_a = zeros(nx,nt);
K_b = zeros(nx,nt);

R_OH = zeros(nx,nt);
R_HA = zeros(nx,nt);
R_Na = zeros(nx,nt);
R_Cl = zeros(nx,nt);
R_C = zeros(nx,nt);
R_H = zeros(nx,nt);
R_A = zeros(nx,nt);

% --- Set IC and BC

G_OH(:,1)= 10000;
J_OH(1,:)= (u_c_ekr(1,:) + u_e_OH(1,:))*10000;

G_HA(:,1)= 10000;
J_HA(1,:)= (u_c_ekr(1,:) + u_e_HA(1,:))*10000;

G_Na(:,1)= 10000;
J_Na(1,:)= (u_c_ekr(1,:) + u_e_Na(1,:))*10000;

G_Cl(:,1)= 10000;
J_Cl(1,:)= (u_c_ekr(1,:) + u_e_Cl(1,:))*10000;

G_C(:,1)= 10000;
J_C(1,:)= (u_c_ekr(1,:) + u_e_C(1,:))*10000;

G_H(:,1)= 10000;
J_H(1,:)= (u_c_ekr(1,:) + u_e_H(1,:))*10000;

G_A(:,1)= 10000;
J_A(1,:)= (u_c_ekr(1,:) + u_e_A(1,:))*10000;

R_OH(:,1) = R_i*(dt)/n;
R_HA(:,1) = R_i*(dt)/n;
R_Na(:,1) = R_i*(dt)/n;
R_Cl(:,1) = R_i*(dt)/n;
R_H(:,1) = R_i*(dt)/n;
R_C(:,1) = R_i*(dt)/n;
R_Cl(:,1) = R_i*(dt)/n;

s_H(:,1) = (z_H^2)*v_H*G_H(:,1);
s_OH(:,1) = (z_OH^2)*v_OH*G_OH(:,1);
s_C(:,1) = (z_C^2)*v_C*G_C(:,1);
Sigma(:,1) = (F^2)*(s_H(:,1) + s_OH(:,1) + s_C(:,1));

for m= 1:nt-1

    G_OH(1,m) =J_OH(1,m);   %--- Upper boundary
    G_HA(1,m) =J_HA(1,m);   %--- Upper boundary
    G_Na(1,m) =J_Na(1,m);   %--- Upper boundary
    G_Cl(1,m) =J_Cl(1,m);   %--- Upper boundary
    G_C(1,m) =J_C(1,m);     %--- Upper boundary
    G_H(1,m) =J_H(1,m);     %--- Upper boundary
    G_A(1,m) =J_A(1,m);     %--- Upper boundary
    
    for i= 2:nx-1
        
        G_OH(i,m+1) = G_OH(i,m) + alpha_OH*(G_OH(i+1,m) -2*G_OH(i,m) + G_OH(i-1,m)) + beta_OH(i,m)*(G_OH(i+1,m) - G_OH(i-1,m)) + R_OH(i,m)/R_D;
        G_HA(i,m+1) = G_HA(i,m) + alpha_HA*(G_HA(i+1,m) -2*G_HA(i,m) + G_HA(i-1,m)) + beta_HA(i,m)*(G_HA(i+1,m) - G_HA(i-1,m)) + R_HA(i,m)/R_D;
        G_Na(i,m+1) = G_Na(i,m) + alpha_Na*(G_Na(i+1,m) -2*G_Na(i,m) + G_Na(i-1,m)) + beta_Na(i,m)*(G_Na(i+1,m) - G_Na(i-1,m)) + R_Na(i,m)/R_D;
        G_Cl(i,m+1) = G_Cl(i,m) + alpha_Cl*(G_Cl(i+1,m) -2*G_Cl(i,m) + G_Cl(i-1,m)) + beta_Cl(i,m)*(G_Cl(i+1,m) - G_Cl(i-1,m)) + R_Cl(i,m)/R_D;        
        G_C(i,m+1) = G_C(i,m) + alpha_C*(G_C(i+1,m) -2*G_C(i,m) + G_C(i-1,m)) + beta_C(i,m)*(G_C(i+1,m) - G_C(i-1,m)) + R_C(i,m)/R_D;
        G_H(i,m+1) = G_H(i,m) + alpha_H*(G_H(i+1,m) -2*G_H(i,m) + G_H(i-1,m)) + beta_H(i,m)*(G_H(i+1,m) - G_H(i-1,m)) + R_H(i,m)/R_D;
        G_A(i,m+1) = G_A(i,m) + alpha_A*(G_A(i+1,m) -2*G_A(i,m) + G_A(i-1,m)) + beta_A(i,m)*(G_A(i+1,m) - G_A(i-1,m)) + R_A(i,m)/R_D;
        

        J_OH(i,m) = (u_c_ekr(i,m) + u_e_OH(i,m))*G_OH(i,m) - alpha_OH*(G_OH(i,m) - G_OH(i-1,m));
        J_HA(i,m) = (u_c_ekr(i,m) + u_e_HA(i,m))*G_HA(i,m) - alpha_HA*(G_HA(i,m) - G_HA(i-1,m));
        J_Na(i,m) = (u_c_ekr(i,m) + u_e_Na(i,m))*G_Na(i,m) - alpha_Na*(G_Na(i,m) - G_Na(i-1,m));
        J_Cl(i,m) = (u_c_ekr(i,m) + u_e_Cl(i,m))*G_Cl(i,m) - alpha_Cl*(G_Cl(i+1,m) - G_Cl(i,m));
        J_C(i,m) = (u_c_ekr(i,m) + u_e_C(i,m))*G_C(i,m) - alpha_C*(G_C(i+1,m) - G_C(i,m));
        J_H(i,m) = (u_c_ekr(i,m) + u_e_H(i,m))*G_H(i,m) - alpha_H*(G_H(i+1,m) - G_H(i,m));
        J_A(i,m) = (u_c_ekr(i,m) + u_e_A(i,m))*G_A(i,m) - alpha_A*(G_A(i+1,m) - G_A(i,m));
        
        
        G_OH(end,m) = J_OH(i,m);    %--- Lower boundary
        G_HA(end,m) = J_HA(i,m);    %--- Lower boundary
        G_Na(end,m) = J_Na(i,m);    %--- Lower boundary
        G_Cl(end,m) = J_Cl(i,m);    %--- Lower boundary
        G_C(end,m) = J_C(i,m);      %--- Lower boundary
        G_H(end,m) = J_H(i,m);      %--- Lower boundary
        G_A(end,m) = J_A(i,m);      %--- Lower boundary
        

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
        K_b(i,m) = (G_Cl(i,m)*G_OH(i,m))/G_Na(i,m);
        R_H(i,m) = (K_H2O(i,m)*G_H(i,m)) + (K_a(i,m)*G_HA(i,m));
        R_OH(i,m) = (K_H2O(i,m)*G_OH(i,m)) + (K_b(i,m)*G_Na(i,m));
        R_Cl(i,m) = (K_b(i,m)*G_Na(i,m));
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