% Define problem parameters
D = 1*10^-9; % Diffusion coefficient
C0 = 0.5; % Initial concentration
L = 0.4; % Length of the domain
T = 40; % Total simulation time
dt = 0.1; % Time step size
dx = 0.01; % Grid spacing

% Create grid
x = 0:dx:L;
N = length(x);

% Initialize concentration
C = C0 * ones(1, N);

% Set up loop for time steps
numSteps = T / dt;
for step = 1:numSteps
    % Calculate concentration gradients
    dC = diff(C) / dx;
    
    % Update concentration using Nernst-Planck equation
    C(2:end-1) = C(2:end-1) + D * dt * (dC(1:end-1) - dC(2:end));
    
    % Apply boundary conditions if necessary
    % For example, C(1) = C0 at left boundary
    
    % Plot concentration at each time step (optional)
    plot(x, C);
    xlabel('Position');
    ylabel('Concentration');
    title(['Time Step: ', num2str(step)]);
    drawnow;
end
