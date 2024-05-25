L = 0.4;                       % length of domain in x direction [m]

I0 = 0.0000829*24*60*60;        % Bq m^-2 day^-1
R1 = (0.693/53.2);              % Decay constant day^-1
tmax = 35;                      % end time [day]

nx = 41;                        % number of nodes in x direction
nt = 50401;                     % number of time steps

dx = L/(nx-1);
dt = tmax/(nt-1);

alpha= 3.5447*10^-9*24*3600;   % Mass advection [m^2/day]
r = alpha*dt/dx^2; rr = 1 - 2*r-R1*dt;

T = 25+273;                            % Absolute temperature [K]
mu_oil = 51*24*3600;              % Oil Viscosity [kg/(m.day)]
mu_water = 0.001*24*3600;         % Water Viscosity [kg/(m.day)]

dEdx = 150;                       % Voltage gradient [V/m]
V = dEdx*L;                       % Voltage [V]
zeta = -0.0027;                   % Zeta potential [V]
u_e = 0;                          % Electromigration mobility [m^2/V.s] 
e_0 = 8.854*10^-12;               % permittivity of free space [F/m]
e_r = 7.5;                        % relative permittivity of clay [F/m]
e_oil = 2.3;                      % Crude oil permittivity [F/m]
e_clay = e_0*e_r;                 % Clay permittivity [F/m]
e_w = 5;                          % Water permittivity [F/m]
n=0.64;                           % Porosity
tau = 0.44;                       % Tortuosity
z = 1;                            % Valency
rho = 1620;                       % Clay Density [kg/m^3]

% Calculated constants
k_eo1 = (e_oil*zeta)*n/(mu_oil);                 % Electroosmotic mobility Alshawabkeh [m^2/V.s]
k_eo2 = (e_oil*zeta)*n/(mu_oil*(tau^2));         % Electroosmotic mobility Vane [m^2/V.s]
k_eo3 = (e_oil*zeta*dEdx)/(mu_oil);              % Electroosmotic mobility Shapiro [S/m]

v0 = 2.5*10^-6; % Convection velocity m day^-1

v1 = u_e + k_eo1;

beta = v*dEdx;

rr1 = beta/(2*dx); rrr = (R1*dt);


J0 = I0/sqrt(R1*alpha); % total inventory of Be-7 in soil

% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);
U = zeros(nx,nt);
S = zeros(nx,nt);
X = zeros(nx,nt);
R = zeros(nx,nt);
E = zeros(nx,nt);
W = zeros(nx,nt);

% --- Set IC and BC
U(:,1)= 2000;
S(:,1)= 2000;
X(:,1)= 4000;
R(:,1)= 6000;
E(:,1)= 8000;
W(:,1)= 10000;

% --- Loop over time steps

for m= 2:nt
    U(1,m) = J0; %--- Upper boundary
    U(end,m) = 0; %--- Lower boundary
    X(1,m) = J0; %--- Upper boundary
    X(end,m) = 0; %--- Lower boundary
    R(1,m) = J0; %--- Upper boundary
    R(end,m) = 0; %--- Lower boundary
    E(1,m) = J0; %--- Upper boundary
    E(end,m) = 0; %--- Lower boundary    
    W(1,m) = J0; %--- Upper boundary
    W(end,m) = 0; %--- Lower boundary 
    
    for i= 2:nx-1
        
        
        U(i,m) = r*U(i-1,m-1)+ rr*U(i,m-1)+ r*U(i+1,m-1);
        
        
    end
end

for m= 2:nt-1

    S(1,m) = I0*J0; %--- Upper boundary
    S(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        
        S(i,m+1) = S(i,m) + r*(S(i+1,m) -2*S(i,m) + S(i-1,m)) - rr1*(S(i+1,m) - S(i-1,m)) - rrr;

    end
end


figure;
plot(t,U(10,:),'-','DisplayName', 'Original FDM');
hold on;

plot(t,S(10,:),'--','DisplayName', 'Optimized FDM');



xlabel('Time');
ylabel('Conc(m)');
set(gca,'YDir','reverse','XAxisLocation','top');
legend();



