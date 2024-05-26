% Geomesh info
L = 0.4;                        % length of domain in x direction [m]       
tmax = 35;                      % end time [day]
nx = 41;                        % number of nodes in x direction
nt = 50401;                     % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);


% Physical info
T = 25 + 273;                     % Temperature [K]
R = 8.314;                        % Gas constant [J/mol.K]
n = 0.64;                         % Porosity
F = 96485;                        % Faraady constant [C/mol]
phi = 25;                         % Electric field [V]
tau = 1.25;                       % Tortuosity
dzdx = 1/tau;                     % divertion field
dEdx = phi/L;                     % Electric field in lenght [V/m]
epsilon = 7*10^10;                % Electrical permitivity [F/m]
mu_oil = 510*24*3600;             % Oil viscosity [kg/(m.day)]
mu_solution = 0.001*24*3600;      % Solution viscosity [kg/(m.day)]
zeta = -0.0027;                   % Zeta potential [V]
K = 0.02;
K_A = 1.75*10^-5;                 % Dissociation constant [mol/m3]
k_i = 0.075;                      % Exprimental constant
coeff = 1/(1+k_i);                % Adsorbing coefficent


% Dimensionless parameters
Pe = 47;
Z = 0.049;
Beta = 967;

% Species Valency
z_H = 1;
z_OH = -1;
z_C = 0;


% Current density
i = F*(z_H+z_OH+z_C);             % there is a flux term at the end as well but it is not calculated since sum of valencies are zero in this phenomena

%-------------------------------------------------------------------------------
% Species Velocities                                                           |
                                                                               |
% Species Mobility                                                             |
v_C = (D_C/(R*T));                % mobility [Hydrocarbon]                     |
v_H = (D_H/(R*T));                % mobility [Hydrogen]                        |
v_OH = (D_OH/(R*T));              % mobility [Hydroxid]                        |
                                                                               |
% Species electromigration velocity                                            |
u_e_H = (v_H*z_H*F*dEdx)/(tau^2);              % electromigration [Hydrogen]   |
u_e_OH = (v_OH*z_OH*F*dEdx)/(tau^2);           % electromigration [Hydroxid]   |
u_e_C = (v_C*z_C*F*dEdx)/(tau^2);              % electromigration [Carbon]     |
                                                                               |
% domain velocity                                                              |
u_x = (epsilon/mu_solution)*(Zeta*dEdx)        % Volumetric Velocity [m3/s]    |
u_c = (1/tau^2)*u_x;                           % Convection Velocity [m3/s]    |
u_eo = ((epsilon*Zeta)/mu_solution)*dEdx;      % Electoosmotic Velocity [m3/s] |
u_s = n*u_c;                                                                   |
%-------------------------------------------------------------------------------

% Counductivity
sigma_surface = 2.74*10^7;        % Conductivity [S/m]


%-------------------------------------------------------------------------------
% Species Diffuision                                                            |
D_H = 3.5447*10^-9*24*3600;       % Mass advection for Hydrogen [m^2/day]       |
D_C = 2.063*10^-9*24*3600;        % Mass advection for Hydrocarbon [m^2/day]    |
D_OH = 0.450*10^-8*24*3600;       % Mass advection for Hydroxid [m^2/day]       |
                                                                                |
% Species Diffuision abberation with coeff                                      |
D_star_H = D_H*coeff*(dt/dx^2);   % Dimensionless of Diffusion                  |
D_star_C = D_C*coeff*(dt/dx^2);   % Dimensionless of Diffusion                  |
D_star_OH = D_OH*coeff*(dt/dx^2); % Dimensionless of Diffusion with K influence |
                                                                                |
% Species Diffuision abberation standalone                                      |
D_prime_H = D_H*(dt/dx^2);        % Dimensionless of Diffusion                  |
D_prime_C = D_C*(dt/dx^2);        % Dimensionless of Diffusion                  |
D_prime_OH = D_OH*(dt/dx^2);      % Dimensionless of Diffusion with K influence |
                                                                                |
% Species Diffusion advection with coeff                                        |
alpha_H = D_star_H/n;             % Diffusion Advection                         |
alpha_C = D_star_C/n;             % Diffusion Advection                         |
alpha_OH = D_star_OH/n;           % Diffusion Advection                         |
                                                                                |
% Species Diffusion advection standalone                                        |
prime_H = D_prime_H/n;            % Diffusion Advection                         |
prime_C = D_prime_C/n;            % Diffusion Advection                         |
prime_OH = D_prime_OH/n;          % Diffusion Advection                         |
%-------------------------------------------------------------------------------

% Double layer
D_l = sqrt((epsilon*R*T)/2*(z^2)*(F^2)*c)


% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);

G_C = zeros(nx,nt);
G_H = zeros(nx,nt);
G_OH = zeros(nx,nt);

Sigma = zeros(nx,nt);
Sigma_ref = ones(nx,nt);
sigma_ref = Sigma_ref*sigma_surface;

s_H = zeros(nx,nt);
s_C = zeros(nx,nt);
s_OH = zeros(nx,nt);

K_H2O = zeros(nx,nt);
K_a = zeros(nx,nt);
K_b = zeros(nx,nt);

% --- Set IC and BC

G_C(:,1)= 10000;
G_C(:,2)= 10000;

G_H(:,1)= 10000;
G_H(:,2)= 10000;

G_OH(:,1)= 2000;
G_OH(:,2)= 2000;


s_H(:,1) = (z_H^2)*v_H*G_H_i;
s_OH(:,1) = (z_OH^2)*v_OH*G_OH_i;
s_C(:,1) = (z_C^2)*v_C*G_C_i;
sigma_coeff = 
Sigma(:,1) = (F^2)*(s_H + s_OH + s_C);

for m= 2:nt-1

    G_C(1,m) =J0; %--- Upper boundary
    G_C(end,m) = 0; %--- Lower boundary

    G_H(1,m) =J0; %--- Upper boundary
    G_H(end,m) = 0; %--- Lower boundary

    G_OH(1,m) =J0; %--- Upper boundary
    G_OH(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        
        G_C(i,m+1) = G_C(i,m) + sub(i,m)*(alpha*(G_C(i+1,m) -2*G_C(i,m) + G_C(i-1,m)) + beta*(G(i+1,m) - G(i-1,m)) + rrr);
        G_H(i,m+1) = G_H(i,m) + sub(i,m)*(alpha_H*(G_H(i+1,m) -2*G_H(i,m) + G_H(i-1,m)) + beta_H*(G_H(i+1,m) - G_H(i-1,m)) + rrr);
        G_OH(i,m+1) = G_OH(i,m) + sub(i,m)*(alpha_OH*(G_OH(i+1,m) -2*G_OH(i,m) + G_OH(i-1,m)) + beta_OH*(G_OH(i+1,m) - G_OH(i-1,m)) + rrr);
        s_H(i,m) = (z_H^2)*v_H*G_H(i,m);
        s_OH(i,m) = (z_OH^2)*v_OH*G_OH(i,m);
        s_C(i,m) = (z_C^2)*v_H*G_C(i,m);
        Sigma(i,m) = (F^2)*(s_C(i,m) + s_H(i,m) + s_OH(i,m)) + Sigma_ref(i,m);
        i_z(i,m) = (-1*Sigma(i,m)*dEdx - F*((z_C*D_C*(G_C(i+1,m) - G_C(i-1,m))) + (z_H*D_H*(G_H(i+1,m) - G_H(i-1,m))) + (z_OH*D_OH*(G_OH(i+1,m) - G_OH(i-1,m)))))/(tau^2);
        if i == 2
            R_prime_H = i_z/F;
            R_H = -1*R_prime_H
        if i == nx-1
            R_prime_OH = i_z/F;
            R_OH = -1*R_prime_OH;
        % K_H2O(i,m) = G_H(i,m)*G_OH(i,m);
        % K_a(i,m) = (G_H(i,m)*G_A(i,m))/G_HA(i,m);
        % K_b(i,m) = (G_B(i,m)*G_OH(i,m))/G_BOH(i,m);
        % R_H = (K_H2O(i,m)*G_H(i,m)) + (K_a(i,m)*G_HA(i,m));
        % R_OH = (K_H2O(i,m)*G_OH(i,m)) + (K_b(i,m)*G_BOH(i,m));
        % R_B = (K_b(i,m)*G_BOH(i,m));
        % R_A = (K_a(i.m)*G_HA(i,m));
        % R_C?
        % 
    end
end