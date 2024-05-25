L = 40;                     % length of domain in x direction [cm]

t_days = 35;                % end time [days]
tmax = t_days*(24*60);      % total time [min]
nx = 41;                    % number of nodes in x direction
nt = tmax +1;               % number of time steps
dx = L/(nx-1);              % position step size [cm]
dt = tmax/(nt-1);           % time step size [min]

% Permitivity equs
F = 96485.33212;            % Faraday constant [C/mol]
z = 1;                      % Valency
T = 25+273;                 % Absolute temperature [K]
e_r = 7.5;                  % relative permittivity of clay [F/m]
e_0 = 8.854*10^-12;         % permittivity of free space [F/m]
e_clay = e_0*e_r;           % Clay permittivity [F/m]
e_oil = 2.2;                % Crude oil permittivity [F/m]
R1 = -200;
c1 = 200;

% constants
n=0.64;                     % Porosity
tau = 0.44;                 % Tortuosity
zeta = -0.0027;             % Zeta potential [V]
dEdx = 1.5;                 % Voltage gradient [V/cm]
I0 = 8000;                  % Initial concentration [mg/Kg]
L1 = (0.693/53.2);          % Same as R let it be for now assumed to be [s^-1]
D0= 3.5447 * 10^-8;         % Diffusion coefficient [m^2.s^-1]
mu_oil = 51;                    % Viscosity [kg/(m.s)]

v = 0;                      % Electromigration velocity [m^2/V.s]

Debye = sqrt((e_oil*T)/(2*(F^2)*(z^2)*I0)); % Debye length                               % reciprocal of the Debye length
k_eo = (e_oil*zeta)*n/(mu_oil);         % Electroosmotic mobility [m^2/V.s] Based on epsilon of both crude oil and clay, might be true form
k_eo2 = (e_w*zeta)*n/(mu_water);

J0 = -D0 - (I0*(v + k_eo)*dEdx);                % total concentration of crude oil in soil
Jw = -D1*c1 - (c1*(u_e + k_eo2)*dEdx);
C0 = 1;
% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);
U = zeros(nx,nt);

% --- Set IC and BC
U(:,1)= 8000;

% --- Loop over time steps

for m= 2:nt
    U(1,m) = I0*Jw; %--- Upper boundary
    U(end,m) = 0; %--- Lower boundary
    
    for i= 2:nx-1
        
        
        U(i,m) = U(i-1,m-1)+ U(i,m-1)+ U(i+1,m-1);
        
        
    end
end

% --- Compare with exact solution at the end of the simulation


plot(t,U(10,:),'--');
xlabel('Time');
ylabel('Conc');
set(gca,'YDir','reverse','XAxisLocation','top');


figure(2)
fig = plot(U(:,1),x');
for k=2:tmax+1
    set(fig,'xdata',U(:,k),'ydata',x')
    set(gca,'YDir','reverse','XAxisLocation','top');
    title('Surface plot of solution.');
    ylabel('Distance (m)');
    pause(.1)
end
