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
K = 0.02;                         % Exprimental Microbal constant
K_A = 1.75*10^-5;                 % Dissociation constant [mol/m3]
k_i = 0.075;                      % Exprimental constant
coeff = 1/(1+k_i);                % Adsorbing coefficent
R_i = (0.693/53.2);               % Initial Reaction flow rate

% Dimensionless parameters
Pe = 47;
Z = 0.049;
Beta = 967;

% Species Valency
z_H = 1;
z_OH = -1;
z_C = 0;
z_HA = 0;
z_BOH = 0;
z_A = -1;
z_B = 1;

% Current density
i = F*(z_H+z_OH+z_C+z_HA+z_BOH+z_A+z_B);
% there is a flux term at the end as well but it is not calculated since sum of valencies are zero in this phenomena


% Counductivity
sigma_surface = 2.74*10^7;        % Conductivity [S/m]


%-------------------------------------------------------------------------------
% Species Diffuision                                                            |
D_H = 3.5447*10^-9*24*3600;       % Mass advection for Hydrogen [m^2/day]       |
D_C = 2.063*10^-9*24*3600;        % Mass advection for Hydrocarbon [m^2/day]    |
D_OH = 0.450*10^-8*24*3600;       % Mass advection for Hydroxid [m^2/day]       |
D_HA = 1.2*10^-8*24*3600;         % Mass advection for Acid [m^2/day]           |
D_A = 1.2*10^-8*24*3600;          % Mass advection for Acid agent [m^2/day]     |
D_BOH = 1.2*10^-8*24*3600;        % Mass advection for Base [m^2/day]           |
D_B = 1.34*10^-8*24*3600;         % Mass advection for Base agent [m^2/day]     |
%                                                                               |
% Species Diffuision abberation with coeff                                      |
D_star_H = D_H*(dt/dx^2);   % Dimensionless of Diffusion                        |
D_star_C = D_C*(dt/dx^2);   % Dimensionless of Diffusion                        |
D_star_OH = D_OH*(dt/dx^2); % Dimensionless of Diffusion                        |
D_star_HA = D_HA*(dt/dx^2); % Dimensionless of Diffusion                        |
D_star_BOH = D_BOH*(dt/dx^2); % Dimensionless of Diffusion                      |
D_star_A = D_A*(dt/dx^2);     % Dimensionless of Diffusion                      |
D_star_B = D_B*(dt/dx^2);     % Dimensionless of Diffuision                     |
%                                                                               |
% Species Diffuision abberation standalone                                      |
D_prime_H = D_H*coeff*(dt/dx^2);        % Dimensionless of Diffusion            |
D_prime_C = D_C*coeff*(dt/dx^2);        % Dimensionless of Diffusion            |
D_prime_OH = D_OH*coeff*(dt/dx^2);      % Dimensionless of Diffusion            |
D_prime_HA = D_HA*coeff*(dt/dx^2);      % Dimensionless of Diffusion            |
D_prime_BOH = D_BOH*coeff*(dt/dx^2);    % Dimensionless of Diffusion            |
D_prime_A = D_A*coeff*(dt/dx^2);        % Dimensionless of Diffusion            |
D_prime_B = D_B*coeff*(dt/dx^2);        % Dimensionless of Diffusion            |
%                                                                               |
% Species Diffusion advection with coeff                                        |
alpha_H = D_star_H/n;             % Diffusion Advection                         |
alpha_C = D_star_C/n;             % Diffusion Advection                         |
alpha_OH = D_star_OH/n;           % Diffusion Advection                         |
alpha_HA = D_star_HA/n;           % Diffusion Advection                         |
alpha_BOH = D_star_BOH/n;         % Diffusion Advection                         |
alpha_A = D_star_A/n;             % Diffusion Advection                         |
alpha_B = D_star_B/n;             % Diffusion Advection                         |
%                                                                               |
% Species Diffusion advection standalone                                        |
alpha_prime_H = D_prime_H/n;            % Diffusion Advection                   |
alpha_prime_C = D_prime_C/n;            % Diffusion Advection                   |
alpha_prime_OH = D_prime_OH/n;          % Diffusion Advection                   |
alpha_prime_HA = D_prime_HA/n;          % Diffusion Advection                   |
alpha_prime_BOH = D_prime_BOH/n;        % Diffusion Advection                   |
alpha_prime_A = D_prime_A/n;            % Diffusion Advection                   |
alpha_prime_B = D_prime_B/n;            % Diffusion Advection                   |
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Species Velocities                                                            |
%                                                                               |
% Species Mobility                                                              |
v_C = (D_C/(R*T));                % mobility [Hydrocarbon]                      |
v_H = (D_H/(R*T));                % mobility [Hydrogen]                         |
v_OH = (D_OH/(R*T));              % mobility [Hydroxid]                         |
v_HA = (D_HA/(R*T));              % mobility [Acid]                             |
v_BOH = (D_BOH/(R*T));            % mobility [Base]                             |
v_A = (D_A/(R*T));                % mobility [A]                                |
v_B = (D_B/(R*T));                % mobility [B]                                |
%                                                                               |
% Species electromigration velocity                                             |
u_e_H = (v_H*z_H*F*dEdx)/(tau^2);              % electromigration [Hydrogen]    |
u_e_OH = (v_OH*z_OH*F*dEdx)/(tau^2);           % electromigration [Hydroxid]    |
u_e_C = (v_C*z_C*F*dEdx)/(tau^2);              % electromigration [Carbon]      |
u_e_HA = (v_HA*z_HA*F*dEdx)/(tau^2);           % electromigration [Acid]        |
u_e_BOH = (v_BOH*z_BOH*F*dEdx)/(tau^2);        % electromigration [Base]        |
u_e_A = (v_A*z_A*F*dEdx)/(tau^2);              % electromigration [A]           |
u_e_B = (v_B*z_B*F*dEdx)/(tau^2);              % electromigration [B]           |
%                                                                               |
% domain velocity                                                               |
u_x = (epsilon/mu_solution)*(zeta*dEdx);       % Volumetric Velocity [m3/s]     |
u_c = ((1/tau^2)*u_x)/(24000000*3600);         % Convection Velocity [m3/s]     |
u_eo = ((epsilon*zeta)/mu_solution)*dEdx;      % Electoosmotic Velocity [m3/s]  |
u_s = n*u_c; %                                                                  |
%                                                                               |
% Species total velocity                                                        |
u_t_C = (u_e_C + u_c)/n; %                                                      |
u_t_H = (u_e_H + u_c)/n; %                                                      |
u_t_OH = (u_e_OH + u_c)/n; %                                                    |
u_t_HA = (u_e_HA + u_c)/n; %                                                    |
u_t_BOH = (u_e_BOH + u_c)/n; %                                                  |
u_t_A = (u_e_A + u_c)/n; %                                                      |
u_t_B = (u_e_B + u_c)/n; %                                                      |
%                                                                               |
% Velocity advection without coefficent                                         |
beta_C = u_t_C*(dt/2*dx); %                                                     |
beta_H = u_t_H*(dt/2*dx); %                                                     |
beta_OH = u_t_OH*(dt/2*dx); %                                                   |
beta_HA = u_t_HA*(dt/2*dx); %                                                   |
beta_BOH = u_t_BOH*(dt/2*dx); %                                                 |
beta_A = u_t_A*(dt/2*dx); %                                                     |
beta_B = u_t_B*(dt/2*dx); %                                                     |
%                                                                               |
% Velocity advection with coefficent abberation                                 |
beta_prime_c = coeff*u_t_C*(dt/2*dx); %                                         |
beta_prime_H = coeff*u_t_H*(dt/2*dx); %                                         |
beta_prime_OH = coeff*u_t_OH*(dt/2*dx); %                                       |
beta_prime_HA = coeff*u_t_HA*(dt/2*dx); %                                       |
beta_prime_BOH = coeff*u_t_BOH*(dt/2*dx); %                                     |
beta_prime_A = coeff*u_t_A*(dt/2*dx); %                                         |
beta_prime_B = coeff*u_t_B*(dt/2*dx); %                                         |
%                                                                               |
%-------------------------------------------------------------------------------

R_D = coeff*(dt)/n;              % Reaction rate Dimensionless factor
% Double layer
% D_l = sqrt((epsilon*R*T)/2*(z^2)*(F^2)*c)
% this formulas exits to confirm there is a dl here, it is not used in anywhere in model

% Flux properties
I0 = 0.0000829*24*60*60;        % Initial concetration
J0 = I0/sqrt(R_i*alpha_C); 
J0 = J0/100;

% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);

sub = zeros(nx,nt);
M_g = exp(-K*t);             % Microbal growth constant 
for w = 1:nx
    sub(w,:) = M_g;
end

G_C = zeros(nx,nt);
G_H = zeros(nx,nt);
G_OH = zeros(nx,nt);
G_HA = zeros(nx,nt);
G_BOH = zeros(nx,nt);
G_A = zeros(nx,nt);
G_B = zeros(nx,nt);

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

R_C = zeros(nx,nt);
R_H = zeros(nx,nt);
R_OH = zeros(nx,nt);
R_HA = zeros(nx,nt);
R_BOH = zeros(nx,nt);
R_A = zeros(nx,nt);
R_B = zeros(nx,nt);

% --- Set IC and BC

G_C(:,1)= 10000;
G_C(:,2)= 10000;

G_H(:,1)= 2000;
G_H(:,2)= 2000;

G_OH(:,1)= 10000;
G_OH(:,2)= 10000;

G_HA(:,1)= 10000;
G_HA(:,2)= 10000;

G_BOH(:,1)= 10000;
G_BOH(:,2)= 10000;

G_A(:,1)= 10000;
G_A(:,2)= 10000;

G_B(:,1)= 10000;
G_B(:,2)= 10000;

R_C(:,1) = R_i*coeff*(dt)/n;
R_OH(:,1) = R_i*coeff*(dt)/n;
R_H(:,1) = R_i*coeff*(dt)/n;
R_HA(:,1) = R_i*coeff*(dt)/n;
R_BOH(:,1) = R_i*coeff*(dt)/n;
R_B(:,1) = R_i*coeff*(dt)/n;
R_A(:,1) = R_i*coeff*(dt)/n;

s_H(:,1) = (z_H^2)*v_H*G_H(:,1);
s_OH(:,1) = (z_OH^2)*v_OH*G_OH(:,1);
s_C(:,1) = (z_C^2)*v_C*G_C(:,1);
Sigma(:,1) = (F^2)*(s_H(:,1) + s_OH(:,1) + s_C(:,1));

for m= 2:nt-1

    G_C(1,m) =J0; %--- Upper boundary
    G_C(end,m) = 0; %--- Lower boundary

    G_H(1,m) =J0; %--- Upper boundary
    G_H(end,m) = 0; %--- Lower boundary

    G_OH(1,m) =J0; %--- Upper boundary
    G_OH(end,m) = 0; %--- Lower boundary

    G_HA(1,m) =J0; %--- Upper boundary
    G_HA(end,m) = 0; %--- Lower boundary

    G_BOH(1,m) =J0; %--- Upper boundary
    G_BOH(end,m) = 0; %--- Lower boundary

    G_A(1,m) =J0; %--- Upper boundary
    G_A(end,m) = 0; %--- Lower boundary

    G_B(1,m) =J0; %--- Upper boundary
    G_B(end,m) = 0; %--- Lower boundary
    for i= 2:nx-1
        
        
        G_C(i,m+1) = G_C(i,m) + (alpha_C*(G_C(i+1,m) -2*G_C(i,m) + G_C(i-1,m)) + beta_C*(G_C(i+1,m) - G_C(i-1,m)) + R_C(i,m)/R_D);
        G_H(i,m+1) = G_H(i,m) + (alpha_H*(G_H(i+1,m) -2*G_H(i,m) + G_H(i-1,m)) + beta_H*(G_H(i+1,m) - G_H(i-1,m)) + R_H(i,m)/R_D);
        G_OH(i,m+1) = G_OH(i,m) + (alpha_OH*(G_OH(i+1,m) -2*G_OH(i,m) + G_OH(i-1,m)) + beta_OH*(G_OH(i+1,m) - G_OH(i-1,m)) + R_OH(i,m)/R_D);
        G_HA(i,m+1) = G_HA(i,m) + (alpha_HA*(G_HA(i+1,m) -2*G_HA(i,m) + G_HA(i-1,m)) + beta_HA*(G_HA(i+1,m) - G_HA(i-1,m)) + R_HA(i,m)/R_D);
        G_BOH(i,m+1) = G_BOH(i,m) + (alpha_BOH*(G_BOH(i+1,m) -2*G_BOH(i,m) + G_BOH(i-1,m)) + beta_BOH*(G_BOH(i+1,m) - G_BOH(i-1,m)) + R_BOH(i,m)/R_D);
        G_A(i,m+1) = G_A(i,m) + (alpha_A*(G_A(i+1,m) -2*G_A(i,m) + G_A(i-1,m)) + beta_A*(G_A(i+1,m) - G_A(i-1,m)) + R_A(i,m)/R_D);
        G_B(i,m+1) = G_B(i,m) + (alpha_B*(G_B(i+1,m) -2*G_B(i,m) + G_B(i-1,m)) + beta_B*(G_B(i+1,m) - G_B(i-1,m)) + R_B(i,m)/R_D);
        
        s_H(i,m) = (z_H^2)*v_H*G_H(i,m);
        s_OH(i,m) = (z_OH^2)*v_OH*G_OH(i,m);
        s_C(i,m) = (z_C^2)*v_H*G_C(i,m);
        Sigma(i,m) = (F^2)*(s_C(i,m) + s_H(i,m) + s_OH(i,m)) + Sigma_ref(i,m);
        
        i_z(i,m) = (-1*Sigma(i,m)*dEdx - F*((z_C*D_C*(G_C(i+1,m) - G_C(i-1,m))) + (z_H*D_H*(G_H(i+1,m) - G_H(i-1,m))) + (z_OH*D_OH*(G_OH(i+1,m) - G_OH(i-1,m)))))/(tau^2);
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
        R_C(i,m) = R_i*coeff*(dt)/n;
        
    end
end


G_C_B = zeros(nx,nt);
G_H_B = zeros(nx,nt);
G_OH_B = zeros(nx,nt);
G_HA_B = zeros(nx,nt);
G_BOH_B = zeros(nx,nt);
G_A_B = zeros(nx,nt);
G_B_B = zeros(nx,nt);

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

for m= 2:nt-1

    G_C_B(1,m) =J0; %--- Upper boundary
    G_C_B(end,m) = 0; %--- Lower boundary

    G_H_B(1,m) =J0; %--- Upper boundary
    G_H_B(end,m) = 0; %--- Lower boundary

    G_OH_B(1,m) =J0; %--- Upper boundary
    G_OH_B(end,m) = 0; %--- Lower boundary

    G_HA_B(1,m) =J0; %--- Upper boundary
    G_HA_B(end,m) = 0; %--- Lower boundary

    G_BOH_B(1,m) =J0; %--- Upper boundary
    G_BOH_B(end,m) = 0; %--- Lower boundary

    G_A_B(1,m) =J0; %--- Upper boundary
    G_A_B(end,m) = 0; %--- Lower boundary

    G_B_B(1,m) =J0; %--- Upper boundary
    G_B_B(end,m) = 0; %--- Lower boundary
    for i= 2:nx-1
        
        
        G_C_B(i,m+1) = G_C_B(i,m) + sub(i,m)*(alpha_prime_C*(G_C_B(i+1,m) -2*G_C_B(i,m) + G_C_B(i-1,m)) + beta_prime_C*(G_C_B(i+1,m) - G_C_B(i-1,m)) + R_C(i,m)/R_D);
        G_H_B(i,m+1) = G_H_B(i,m) + sub(i,m)*(alpha_prime_H*(G_H_B(i+1,m) -2*G_H_B(i,m) + G_H_B(i-1,m)) + beta_prime_H*(G_H_B(i+1,m) - G_H_B(i-1,m)) + R_H(i,m)/R_D);
        G_OH_B(i,m+1) = G_OH_B(i,m) + sub(i,m)*(alpha_prime_OH*(G_OH_B(i+1,m) -2*G_OH_B(i,m) + G_OH_B(i-1,m)) + beta_prime_OH*(G_OH_B(i+1,m) - G_OH_B(i-1,m)) + R_OH(i,m)/R_D);
        G_HA_B(i,m+1) = G_HA_B(i,m) + sub(i,m)*(alpha_prime_HA*(G_HA_B(i+1,m) -2*G_HA_B(i,m) + G_HA_B(i-1,m)) + beta_prime_HA*(G_HA_B(i+1,m) - G_HA_B(i-1,m)) + R_HA(i,m)/R_D);
        G_BOH_B(i,m+1) = G_BOH_B(i,m) + sub(i,m)*(alpha_prime_BOH*(G_BOH_B(i+1,m) -2*G_BOH_B(i,m) + G_BOH_B(i-1,m)) + beta_prime_BOH*(G_BOH_B(i+1,m) - G_BOH_B(i-1,m)) + R_BOH(i,m)/R_D);
        G_A_B(i,m+1) = G_A_B(i,m) + sub(i,m)*(alpha_prime_A*(G_A_B(i+1,m) -2*G_A_B(i,m) + G_A_B(i-1,m)) + beta_prime_A*(G_A_B(i+1,m) - G_A_B(i-1,m)) + R_A(i,m)/R_D);
        G_B_B(i,m+1) = G_B_B(i,m) + sub(i,m)*(alpha_prime_B*(G_B_B(i+1,m) -2*G_B_B(i,m) + G_B_B(i-1,m)) + beta_prime_B*(G_B_B(i+1,m) - G_B_B(i-1,m)) + R_B(i,m)/R_D);
        
        s_H(i,m) = (z_H^2)*v_H*G_H_B(i,m);
        s_OH(i,m) = (z_OH^2)*v_OH*G_OH_B(i,m);
        s_C(i,m) = (z_C^2)*v_H*G_C_B(i,m);
        Sigma(i,m) = (F^2)*(s_C(i,m) + s_H(i,m) + s_OH(i,m)) + Sigma_ref(i,m);
        
        i_z(i,m) = (-1*Sigma(i,m)*dEdx - F*((z_C*D_C*(G_C_B(i+1,m) - G_C(i-1,m))) + (z_H*D_H*(G_H_B(i+1,m) - G_H_B(i-1,m))) + (z_OH*D_OH*(G_OH_B(i+1,m) - G_OH_B(i-1,m)))))/(tau^2);
        if i == 2
            R_prime_H = i_z(i,m)/F;
            R_H(i,m) = -1*R_prime_H;
        end
        if i == nx-1
            R_prime_OH = i_z(i,m)/F;
            R_OH(i,m) = -1*R_prime_OH;
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
pH_scale = linspace(1,40,41);
xl = [0,5,10,15,20,25,30,35];
yl = [10000,7900,7100,6000,5700,5500,5400,5100];
figure;

hold on;

plot(t,G_C(10,:),'--','DisplayName', 'Hydrocarbon');

% plot(t,G_H(10,:),'-','DisplayName', 'Hydrogen');

plot(t,G_OH(10,:),'-','DisplayName', 'Hydroxid');

plot(t,G_BOH(10,:),'-','DisplayName', 'Base');

plot(t,G_HA(10,:),'-','DisplayName', 'Acid');

scatter(xl,yl, 'DisplayName', 'Expriment Data');

xlabel('Time');
ylabel('Conc(mg/kg)');

legend();


figure;
hold on;
plot(pH_scale,pH(:,50400),'--','DisplayName', 'pH')

xlabel('Length (cm)');
ylabel('pH level');

legend();
