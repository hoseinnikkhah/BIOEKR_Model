% Define function to solve transport equation
function [c, E] = solve_transport(c0, D1, D2, F, zi, zy, Di, Dy, sigma, L, Nx, dx)
    % Initialize concentration and electric field arrays
    c = c0;
    E = zeros(Nx, 1);
    
    % Iterate through time steps
    for t = 1:1000
        % Calculate electric field
        E(1:end-1) = -sigma*D1*c(1:end-1,:)*zi'/F;
        E(end) = 0; % Boundary condition at x=L
        
        % Calculate fluxes
        J = zeros(Nx, Ni+Ny);
        for i = 1:Ni
            J(:,i) = -F*zi(i)*Di(i)*D1*c(:,i);
        end
        for y = 1:Ny
            J(:,Ni+y) = -F*zy(y)*Dy(y)*D1*c(:,Ni+y);
        end
        
        % Update concentration
        c = c - D2*J*dx + D1*diag(E)*c*dx;
    end
    % Solve transport equation
    [c, E] = solve_transport(c0, D1, D2, F, zi, zy, Di, Dy, sigma, L, Nx, dx);
    % Plot results
    x = linspace(0, L, Nx);
    figure;
    plot(x, c(:,1), 'b-', x, c(:,2), 'r-', x, c(:,3), 'g-', x, E, 'k-');
    xlabel('x');
    ylabel('Concentration/Electric field');
    legend('Species 1', 'Species 2', 'Species y', 'Electric field');
end
