L = 40;                           % length of domain in x direction [cm]
t_obs = 35;                       % total time [days]
tmax = t_obs*(24*60);             % total time [min]
nx = 41;                          % number of nodes in x direction
nt = tmax +1;                     % number of time steps
dx = L/(nx-1);                    % position step size [cm]
dt = tmax/(nt-1);                 % time step size [min]

I0 = 6000;
n = 0.64;
R1 = (0.693);
D0 = 3.5447*10^-9;
D1 = 9.3100*10^-8;

u_e = 0;                          % Electromigration mobility [m^2/V.s] 
dEdx = 1.5;                       % Voltage gradient [V/cm]
zeta = -0.0027;                   % Zeta potential [V]
mu_oil = 0.51;                      % Oil Viscosity [kg/(cm.s)]
e_oil = 0.23;                      % Crude oil permittivity [F/m]

k_eo3 = (e_oil*zeta)*dEdx/(mu_oil);              % Electroosmotic mobility Shapiro [m^2/V.s]

alpha= D0*10000;
r = alpha*dt/(n*dx^2); rr = ((u_e + k_eo3)*dEdx*dt)/(2*n*dx); rrr = (R1*dt)/n;

J0 = -(D1*I0) - I0*(u_e + k_eo3)*dEdx; % total inventory of Be-7 in soil

v = 2.5*10^-6; % Convection velocity m day^-1
v0 = (v)*(dt/dx);
% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);

U = zeros(nx,nt);

% --- Set IC and BC
U(:,1)= 6000;


for m= 2:nt-1
    U(1,m) = 0; %--- Upper boundary
    U(end,m) = 1*J0; %--- Lower boundary


    for i= 2:nx-1
        
        
        U(i,m+1) = U(i,m) + r*(U(i+1,m) -2*U(i,m) + U(i-1,m)) + rr*(U(i+1,m) - U(i-1,m)) + rrr;
        % original cannot resolve correctly
        

    end
end


figure;

% Plot U
plot(t, U(20,:), '-', 'DisplayName', 'U at 20');
hold on;

% Plot S
plot(t, U(18,:), '-', 'DisplayName', 'S at 18');

% Plot K
plot(t, U(10,:), '-', 'DisplayName', 'K at 30');

hold off;

xlabel('Time');
ylabel('Conc(M)');
set(gca, 'YDir', 'reverse', 'XAxisLocation', 'top');
title('Comparison of U, S, and K');
legend();

