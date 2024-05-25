% Define grid
N = 100; % number of nodes
L = 87; % domain length
dx = L/N; % grid spacing
x = linspace(0, L, N+1); % x grid points

% Define boundary conditions
Ci_left = 0; % Dirichlet boundary condition at left boundary
Ci_right = 0.1120; % Dirichlet boundary condition at right boundary
Ci_top = 0; % Dirichlet boundary condition at top boundary
Ci_bottom = 0; % Dirichlet boundary condition at bottom boundary

% Define parameters
Di = 1; % effective diffusion coefficient
ui = 1; % effective ionic mobility coefficient
ke = 1; % electroosmotic permeability coefficient
ce = 1; % concentration of electrolyte
kh = 1; % hydraulic conductivity coefficient
E = 20; % electric potential
DeltaH = 0; % hydraulic head

% Discretize the equation
A = zeros(N-1); % system matrix
B = zeros(N-1, 1); % right-hand side vector
for i = 2:N-2
    A(i, i-1) = -Di/dx^2 - ui*x(i)/(2*dx);
    A(i, i) = 2*Di/dx^2 + ke*ce*E + kh*DeltaH;
    A(i, i+1) = -Di/dx^2 + ui*x(i)/(2*dx);
    B(i) = 0; % right-hand side values
end

% Set boundary conditions
A(1, 1) = 2*Di/dx^2 + ke*ce*E + kh*DeltaH;
A(1, 2) = -2*Di/dx^2 + ui*x(1)/(2*dx);
B(1) = Ci_left;

A(N-1, N-1) = 2*Di/dx^2 + ke*ce*E + kh*DeltaH;
A(N-1, N-2) = -2*Di/dx^2 - ui*x(N)/(2*dx);
B(N-1) = Ci_right;

% Solve the system of equations
Ci = [Ci_left; A\B; Ci_right];
Ci = [Ci_top, Ci(2:end-1)', Ci_bottom]; % add boundary values to Ci

% Plot the solution
figure;
plot(x, Ci);
xlabel('x');
ylabel('Ci');
