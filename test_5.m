% Parameters
L = 40;        % Length of the spatial domain
T = 3000;        % Total simulation time

nx = 41;      % Number of spatial points
nt = 3001;     % Number of time steps

dx = L / (nx-1);    % Spatial step size
dt = T / (nt-1);    % Time step size

K = 1;        % Some constant (adjust as needed)
u = 1;        % convection velocity
k = 1;        % Some constant (adjust as needed)
R = 1;

x = linspace(0,L,nx);
t = linspace(0,tmax,nt);
U = zeros(nx,nt);


% Initial condition
c = zeros(nx, 1);
c0 = 2000;       % Initial concentration
c(1) = c0;

% FTCS method
for m= 2:nt
    U(1,m) = J0; %--- Upper boundary
    U(end,m) = 0; %--- Lower boundary
    
    for i= 2:nx-1
        
        
        U(i,m) = r*U(i-1,m-1)+ rr*U(i,m-1)+ r*U(i+1,m-1);
        
        
    end
end

% Plot the result
x = linspace(0, L, N);
plot(x, c);
xlabel('Position (x)');
ylabel('Concentration (c)');
title('FTCS Solution for Mass Transfer PDE');
