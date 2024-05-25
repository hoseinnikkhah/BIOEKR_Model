% Geomesh info
L = 0.3;                          % length of domain in x direction [m]       
tmax = 150;                       % end time [day]
nx = 31;                          % number of nodes in x direction [cm]
nt = 216001;                      % number of time steps [min]
dx = L/(nx-1);
dt = tmax/(nt-1);

I0 = 0.0000829*24*60*60;
R1 = (0.693/53.2);


% Physical info
T_i = 15 + 273;                     % Diffusion temp [Kelvin]
T_D = 22 + 273;                     % Working temp [Kelvin]
tau = 0.6;                          % Tortuosity
n = 0.25;                           % Porosity
V_Cathod = 0;                       % Electric potential at the Cathode [V]
V_anode = 110;                      % Electric potential at the Anode [V]
R_i = 5*10^-7;                      % dissolution rate coefficient [1/s]
sigma_surface = 0.0013;             % Conductivity
F = 96485;                          % Faraady constant [C/mol]
epsilon_water = 80;                 % dielectric constant [water]
mu_water = (1.002/60)*10^-3;        % Water Viscosity [N⋅s/m^2]
zeta = -0.0027;                     % Zeta potential [V]
R = 8.314;                          % Gas constant [J/mol.K]
e = 43;                             % dielectric constant of the porous medium at 0.6 vulumetric water content
dEdx = (V_anode-V_Cathod)/L;
% Electroosmotic permeability
k_eo = ((n*epsilon_water)/mu_water)*zeta;

% ion charge and Diffusion of solvents in order
z_i = [-1,0,0,0,0,0,1,-1,1,-2,-1,-1,0,2,1,0];
D_i = [8.36*10^-10,6.50*10^-10,7.35*10^-10,8.63*10^-10,1.09*10^-9,1.67*10^-9,7.53*10^-9,1.64*10^-9,1.08*10^-9,7.73*10^-10,9.55*10^-10,9.72*10^-10,5.45*10^-10,6.41*10^-10,4.10*10^-10,3.61*10^-10];
% Lactate,PCE,TCE,DCE,VC,Ethene,H+,Cl-,Na+,CO3,HCO3,NaCO3,NaHCO3,Ca,CaHCO3,CaCO3

D_i_star = (D_i./tau).*(dt/dx^2);
D_i_star1 = (D_i./tau).*(dt/(2*dx));


J0 = I0/sqrt(R1*D_i_star(1,1)); 
J0 = J0/100;



% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);
G = zeros(nx,nt);

% --- Set IC and BC

G(:,1)= 10000;
G(:,2)= 10000;

for m= 2:nt-1

    G(1,m) =J0; %--- Upper boundary
    G(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        
        G(i,m+1) = G(i,m) + D_i_star(1,1)*(G(i+1,m) -2*G(i,m) + G(i-1,m)) + n*(z_i(1,1)/(R*T_i))*D_i_star1(1,1)*dEdx*(G(i+1,m) - G(i-1,m)) - n*k_eo*(dt/(2*dx))*(G(i+1,m) - G(i-1,m)) + R_i;
        
    end
end

plot(t,G(10,:),'-','DisplayName', 'New FDM for Hydrocarbon');