% One-Dimensional Finite Element Method

% Define the problem parameters
L = 1; % Length of the domain
N = 100; % Number of elements
h = L/N; % Element size

% Define the basis functions (shape functions)
phi = @(x) [1-x/h, x/h];

% Define the element stiffness matrix
k = @(x) [1 -1; -1 1];

% Define the element load vector
f = @(x) [x; x^2]; % Example load function

% Assemble the global stiffness matrix and load vector
K = zeros(N+1);
F = zeros(N+1, 1);

for i = 1:N
    x = (i-1)*h; % Element starting point
    Ke = k(x+h/2)*h; % Element stiffness matrix
    Fe = f(x+h/2)*h; % Element load vector
    
    % Assemble the element contributions into the global system
    indices = i:i+1;
    K(indices, indices) = K(indices, indices) + Ke;
    F(indices) = F(indices) + Fe;
end

% Apply boundary conditions
% Modify the appropriate rows and columns of K and F

% Apply Dirichlet boundary conditions
K(1, :) = 0;
K(1, 1) = 1;
F(1) = 0;

% Apply Neumann boundary conditions
% Modify the appropriate entries in F

% Apply initial conditions
U0 = zeros(N+1, 1); % Initial concentration values
U = U0;

% Solve the linear system at each time step
T = 1; % Total simulation time
dt = 0.01; % Time step size
numSteps = T/dt; % Number of time steps

for step = 1:numSteps
    % Solve the linear system
    U = K\F;
    
    % Update the load vector for the next time step
    F = zeros(N+1, 1);
    
    % Apply Neumann boundary conditions
    % Modify the appropriate entries in F
    
    % Apply initial conditions for the next time step
    U(1) = U0(1); % Assuming the initial condition is known at the left boundary
    
    % Plot the solution at each time step (optional)
    x = linspace(0, L, N+1);
    plot(x, U, 'b-');
    xlabel('x');
    ylabel('u(x)');
    title(['Finite Element Solution at Time t = ', num2str(step*dt)]);
    drawnow;
end
