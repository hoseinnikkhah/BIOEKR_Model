
% Geometry parameters
L = 40;                       % length of domain in x direction [cm]
t_obs = 35;                   % total time [days]
tmax = t_obs*(24);         % total time [min]
nx = 41;                      % number of nodes in x direction
nt = tmax +1;                 % number of time steps
dx = L/(nx-1);                % position step size [cm]
dt = tmax/(nt-1);             % time step size [min]

% Global constants
R = 8.314;                     % Gas constant [J/mol.K]
F = 96485.33212;               % Faraday constant [C/mol]

% Physical parameters (are given)
D0 = 3.5447*10^-9;                       % Diffusion (Crud on in clay) [m^2/s] 
D1 = 9.3100*10^-9;                       % Diffusion (H+ on in clay) [m^2/s]
r0 = D0*((1/(dx^2))*10000)*(dt/60);      % Crude oil Diffusion [dimensionless]
r1 = D1*((1/(dx^2))*10000)*(dt/60);      % H+ Diffusion [dimensionless]
rr = 1 -2*r0 - (L1*dt);                 % there are no L1 so this is different from ref code
rr0 = (((u_e+k_eo3)*dEdx)/n)*(dt/2*dx);
v0 = (((u_e+k_eo3)*dEdx)/n)*(dt/dx);

% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);

U = zeros(nx,nt);
K = zeros(nx,nt);
S = zeros(nx,nt);

% --- Set IC and BC
U(:,1) = 0;
S(:,1) = 0;
K(:,1) = 0;

% --- Loop over time steps

for m= 2:nt
    U(1,m) = J0; %--- Upper boundary
    U(end,m) = 0; %--- Lower boundary

    S(1,m) = J0; %--- Upper boundary
    S(end,m) = 0; %--- Lower boundary

    K(1,m) = J0; %--- Upper boundary
    K(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        
        U(i,m) = r0*U(i-1,m-1) + rr*U(i,m-1) + r0*U(i+1,m-1);
        % original cannot resolve correctly
        S(i,m) = r0 * S(i-1,m-1) + rr * S(i,m-1) + r0 * S(i+1,m-1) + (u_e + k_eo3) * dEdx * (S(i+1,m-1) - S(i-1,m-1)) - L1 * S(i,m-1);
        %S(i,m) = r0*S(i-1,m-1) + rr*S(i,m-1) + r0*U(i+1,m-1) + v0*(S(i-1,m-1) - S(i,m-1));
        %based on above with v term included
        %K(i,m) = -K(i,m+1) - (r0)*(K(i+1,m) - 2*K(i,m) + K(i-1,m)) + rr0*(K(i+1,m)-K(i+1,m));
        % new formula based on summery
        K(i,m) = r0*K(i-1,m-1) + rr0*K(i,m-1) + r0*K(i+1,m-1) - L1*dt*K(i,m-1);
        

        
    end
end

% Plotting in same figure to compare
figure;

% Plot U
plot(t, U(1,:), '--', 'DisplayName', 'U');
hold on;

% Plot S
plot(t, S(5,:), '--', 'DisplayName', 'S');

% Plot K
plot(t, K(10,:), '--', 'DisplayName', 'K');

hold off;

xlabel('Time');
ylabel('Conc(M)');
set(gca, 'YDir', 'reverse', 'XAxisLocation', 'top');
title('Comparison of U, S, and K');
legend();
