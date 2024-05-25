% Geomesh info
L = 0.335;                      % length of domain in x direction [m]       
tmax = 35;                      % end time [day]
nx = 35;                        % number of nodes in x direction
nt = 50401;                     % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);

I0 = 0.0000829*24*60*60;        % Initial concetration
R1 = (0.693/53.2);              % Reaction flow rate

% Reaction rate is not constant all the time and is different for each
% species, in order to calculate this factor we need to first calculate
% several constants or factors and then apply it in adequte formula
% followin code is dedicated to this part and might not apply to whole
% boundary condition

% Physical info
n = 0.41;                         % Porosity
tau = 2.56;                       % Tortuosity
z1 = 1;                           % Valancy constant (H+)
z2 = 0;                           % valancy constant (Hydrocarbon)
z3 = -1;                          % valancy constant (OH-)
sigma_surface = 2.74*10^7;        % Conductivity [S/m]
k_i = 0.075;
coeff = 1/(1+k_i);                % Adsorbing coefficent
R = 8.314;                        % Gas constant [J/mol.K]



% H+ Diffusion
D_H = 3.5447*10^-9*24*3600;       % Mass advection [m^2/day]
D_star_H = D_H*coeff*(dt/dx^2);   % Dimensionless of Diffusion
D_star1_H = D_H*(dt/dx^2);
alpha_H = D_star_H/n;             % Diffusion Advection 
alpha1_H = D_star1_H/n;

% Hydrocarbon Diffusion
D = 2.063*10^-9*24*3600;          % Mass advection [m^2/day]
D_star = D*coeff*(dt/dx^2);       % Dimensionless of Diffusion
D_star1 = D*(dt/dx^2);
alpha = D_star/n;                 % Diffusion Advection 
alpha1 = D_star1/n;

% OH- Diffusion
D_OH = 0.450*10^-8*24*3600;       % Mass advection [m^2/day]
D_star_OH = D_OH*coeff*(dt/dx^2); % Dimensionless of Diffusion
D_star1_OH = D_OH*(dt/dx^2);
alpha_OH = D_star_OH/n;           % Diffusion Advection 
alpha1_OH = D_star1_OH/n;


% Inputs
dEdx = 150;                       % Eletric field [V/m]
epsilon = 2.1;                    % dielectric constant [hydrocarbon]
epsilon_H = 1.2;                  % dielectric constant [Hydrogen]
epsilon_water = 80;               % dielectric constant [water]
epsilon_OH = 2.21;                % dielectric constant [OH]
zeta = -0.0027;                   % Zeta potential [V]
T = 25+273;                       % Absolute temperature [K]
mu_water = 0.001*24*3600;         % Water Viscosity [kg/(m.day)]
mu_oil = 510*24*3600;             % Oil viscosity [kg/(m.day)]
F = 96485;                        % Faraady constant [C/mol]
K = 0.02;

% Electroosmotic Values
u_star = ((epsilon*zeta)/mu_oil)*(dEdx/(tau^2));                 % Hydrocarbon
u_star_H = ((epsilon_H*zeta)/mu_water)*(dEdx/(tau^2));           % H+
u_star_water = ((epsilon_water*zeta)/mu_water)*(dEdx/(tau^2));   % Water
u_star_OH = ((epsilon_OH*zeta)/mu_water)*(dEdx/(tau^2));         % OH-

% Electromigration values for Hydrocarbon
v = (D/(R*T));                            % mobility [Hydrocarbon]
u_e = (v*z2*F*dEdx)/(tau^2);              % electromigration [Hydrocarbon]
u_total = coeff*(u_e - u_star);           % Total mobility advection [Hydrocarbon]
beta = (u_total/n)*(dt/2*dx);             % Dimensionless mobility advection [Hydrocarbon]

% Electromigration values for H+
v_H = (D_H/(R*T));                        % mobility [H+]
u_e_H = (v_H*z1*F*dEdx)/(tau^2);          % electromigration [H+]
u_total_H = coeff*(u_e_H - u_star_H);     % Total mobility advection [H+]
beta_H = (u_total_H/n)*(dt/2*dx);         % Dimensionless mobility advection [H+]

% Electromigration values for H+
v_OH = (D_OH/(R*T));                      % mobility [OH-]
u_e_OH = (v_OH*z1*F*dEdx)/(tau^2);        % electromigration [OH-]
u_total_OH = coeff*(u_e_OH - u_star_OH);  % Total mobility advection [OH-]
beta_OH = (u_total_OH/n)*(dt/2*dx);       % Dimensionless mobility advection [OH-]


rrr = coeff*(R1*dt)/n;
u_total1_H = (u_e_H - u_star_H);
beta1_H = (u_total1_H/n)*(dt/2*dx);

J0 = I0/sqrt(R1*alpha); 
J0 = J0/100;


% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);
sub2 = exp(-K*t);
sub = zeros(nx,nt);

for ll = 1:nx
    sub(ll,:) = sub2;
end

ix_H = zeros(nx,nt);
ix_OH = zeros(nx,nt);
Sigma_ref = sigma_surface * ones(nx,nt);
Sigma = zeros(nx,nt);

G = zeros(nx,nt);
G1 = zeros(nx,nt);
G_H = zeros(nx,nt);
G_OH = zeros(nx,nt);
B = zeros(nx,nt);

% --- Set IC and BC

G(:,1)= 10000;
G(:,2)= 10000;

G_H(:,1)= 10000;
G_H(:,2)= 10000;

G_OH(:,1)= 2000;
G_OH(:,2)= 2000;

B(:,1)= 10000;
B(:,2)= 10000;

% --- Loop over time steps

for m= 2:nt-1

    G(1,m) =J0; %--- Upper boundary
    G(end,m) = 0; %--- Lower boundary

    %G_H(1,m) =J0; %--- Upper boundary
    %G_H(end,m) = 0; %--- Lower boundary

    G_OH(1,m) =J0; %--- Upper boundary
    G_OH(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        
        G(i,m+1) = G(i,m) + sub(i,m)*(alpha*(G(i+1,m) -2*G(i,m) + G(i-1,m)) + beta*(G(i+1,m) - G(i-1,m)) + rrr);
        G_H(i,m+1) = G_H(i,m) + sub(i,m)*(alpha_H*(G_H(i+1,m) -2*G_H(i,m) + G_H(i-1,m)) + beta_H*(G_H(i+1,m) - G_H(i-1,m)) + rrr);
        G_OH(i,m+1) = G_OH(i,m) + sub(i,m)*(alpha_OH*(G_OH(i+1,m) -2*G_OH(i,m) + G_OH(i-1,m)) + beta_OH*(G_OH(i+1,m) - G_OH(i-1,m)) + rrr);
        Sigma(i,m) = (F^2)*(z1^2)*v*G(i,m) + Sigma_ref(i,m);
        ix_H(i,m) = (-1*Sigma(i,m) - F*z1*D*(G(i+1,m) - G(i-1,m)))/(tau^2);
        ix_OH(i,m) = (-1*Sigma(i,m) - F*z3*D*(G(i+1,m) - G(i-1,m)))/(tau^2);
        
    end
end

R_H = ix_H/F;
R_OH = ix_OH/F;

J_H = u_total_H*I0 + R_H;
J_OH = u_total_H*I0 + R_OH;
rrr1 = coeff*(R_H*dt)/n;


G1(:,1)= 10000;
G1(:,2)= 10000;

G_H1(:,1)= 10000;
G_H1(:,2)= 10000;

G_OH1(:,1)= 2000;
G_OH1(:,2)= 2000;

B1(:,1)= 10000;
B1(:,2)= 10000;

for m= 2:nt-1

    G1(1,m) =J0; %--- Upper boundary
    G1(end,m) = 0; %--- Lower boundary

    %G_H(1,m) =J0; %--- Upper boundary
    %G_H(end,m) = 0; %--- Lower boundary

    %G_OH1(1,m) =J0; %--- Upper boundary
    %G_OH1(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        
        G1(i,m+1) = G1(i,m) + sub(i,m)*(alpha*(G1(i+1,m) -2*G1(i,m) + G1(i-1,m)) + beta*(G1(i+1,m) - G1(i-1,m)) + rrr1(i,m));
        %G_H1(i,m+1) = G_H1(i,m) + sub(i,m)*(alpha_H*(G_H(i+1,m) -2*G_H1(i,m) + G_H1(i-1,m)) + beta_H*(G_H(i+1,m) - G_H1(i-1,m)) + rrr1);
        %G_OH1(i,m+1) = G_OH1(i,m) + sub(i,m)*(alpha_OH*(G_OH1(i+1,m) -2*G_OH1(i,m) + G_OH1(i-1,m)) + beta_OH*(G_OH1(i+1,m) - G_OH1(i-1,m)) + rrr1);
        %Sigma(i,m) = (F^2)*(z1^2)*v*G(i,m) + Sigma_ref(i,m);
        %ix_H(i,m) = (-1*Sigma(i,m) - F*z1*D*(G(i+1,m) - G(i-1,m)))/(tau^2);
        %ix_OH(i,m) = (-1*Sigma(i,m) - F*z3*D*(G(i+1,m) - G(i-1,m)))/(tau^2);
        
    end
end


for m= 2:nt-1

    G(1,m) =J0; %--- Upper boundary
    G(end,m) = 0; %--- Lower boundary

    G_H(1,m) =J_H(2,50400); %--- Upper boundary
    G_H(end,m) = 0; %--- Lower boundary

    G_OH(1,m) =J0; %--- Upper boundary
    G_OH(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        
        G(i,m+1) = G(i,m) + sub(i,m)*(alpha*(G(i+1,m) -2*G(i,m) + G(i-1,m)) + beta*(G(i+1,m) - G(i-1,m)) + rrr);
        G_H(i,m+1) = G_H(i,m) + sub(i,m)*(alpha_H*(G_H(i+1,m) -2*G_H(i,m) + G_H(i-1,m)) + beta_H*(G_H(i+1,m) - G_H(i-1,m)) + rrr);
        G_OH(i,m+1) = G_OH(i,m) + sub(i,m)*(alpha_OH*(G_OH(i+1,m) -2*G_OH(i,m) + G_OH(i-1,m)) + beta_OH*(G_OH(i+1,m) - G_OH(i-1,m)) + rrr);
        Sigma(i,m) = (F^2)*(z1^2)*v*G(i,m) + Sigma_ref(i,m);
        ix_H(i,m) = (-1*Sigma(i,m) - F*z1*D*(G(i+1,m) - G(i-1,m)))/(tau^2);
        ix_OH(i,m) = (-1*Sigma(i,m) - F*z3*D*(G(i+1,m) - G(i-1,m)))/(tau^2);
        
    end
end



for m= 2:nt-1

    B(1,m) =J0; %--- Upper boundary
    B(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        B(i,m+1) = B(i,m) + alpha1_H*(B(i+1,m) -2*B(i,m) + B(i-1,m)) + beta1_H*(B(i+1,m) - B(i-1,m));

        
    end
end
xl = [0,5,10,15,20,25,30,35];
yl = [10000,7900,7100,6000,5700,5500,5400,5100];
figure;

hold on;

plot(t,B(10,:),'--','DisplayName', 'Original FDM for Hydrocarbon');

plot(t,G(10,:),'-','DisplayName', 'New FDM for Hydrocarbon');

plot(t,G1(10,:),'-','DisplayName', 'New FDM for Hydrocarbon (modifed R');

%plot(t,G_OH(10,:),'-','DisplayName', 'New FDM for OH');

%plot(t,G_H(10,:),'-','DisplayName', 'New FDM for H');

scatter(xl,yl, 'DisplayName', 'Expriment Data');

xlabel('Time');
ylabel('Conc(mg/kg)');

legend();
