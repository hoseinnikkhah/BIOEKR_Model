% Define parameters
L = 1; % Length of domain
Nx = 100; % Number of grid points
dx = L/Nx; % Grid spacing
F = 1; % Faraday's constant
Ni = 2; % Number of species i
Ny = 1; % Number of species y
zi = [-1, 1]; % Valence of species i
zy = [1]; % Valence of species y
Di = [1e-5, 1e-6]; % Diffusion coefficient of species i
Dy = [1e-4]; % Diffusion coefficient of species y
sigma = 1; % Conductivity
c0 = zeros(Nx, Ni+Ny); % Initial concentration
c0(:,1) = 1; % Initial concentration of species 1

% Define finite difference matrices
D2 = gallery('tridiag', Nx, 1, -2, 1)/(dx^2); % Second-order central difference matrix
D1 = gallery('tridiag', Nx, -1, 0, 1)/(2*dx); % First-order central difference matrix
D1(1,:) = 0; % Boundary condition at x=0
D1(end,:) = 0; % Boundary condition at x=L

