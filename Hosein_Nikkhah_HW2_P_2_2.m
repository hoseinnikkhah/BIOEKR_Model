function Hosein_Nikkhah_HW2_P_2_2()
    % Generate a random 10x10 matrix
    A = rand(10);

    % Initial guess for eigenvector and eigenvalue
    v_guess = rand(10, 1);
    lambda_guess = rand();

    % Tolerance for convergence
    tolerance = 1e-8;

    % Perform the Newton-Raphson iteration
    [v, lambda] = newton_raphson_iteration(v_guess, lambda_guess, A, tolerance);

    % Display the results
    disp('Eigenvalue:');
    disp(lambda);
    disp('Eigenvector:');
    disp(v);
    % Check against MATLAB's eig function
    eigenvalue = lambda;   
    [matlab_eigenvectors, matlab_eigenvalues] = eig(A);
    [~, index] = min(abs(diag(matlab_eigenvalues) - eigenvalue));
    matlab_eigenvector = matlab_eigenvectors(:, index);

    disp('MATLAB eig Function Result:');
    disp('Eigenvalue:');
    disp(matlab_eigenvalues(index, index));
    disp('Eigenvector:');
    disp(matlab_eigenvector);      
end

function [v, lambda] = newton_raphson_iteration(v_guess, lambda_guess, A, tolerance)
    v = v_guess;
    lambda = lambda_guess;

    % Maximum number of iterations (adjust as needed)
    max_iterations = 100;

    for iteration = 1:max_iterations
        % Calculate the Jacobian matrix
        J = jacobian_matrix(v, lambda, A);

        % Calculate the function values
        F = Hosein_Nikkhah_HW2_P2_1_1(v, lambda, A);

        % Update using the Newton-Raphson formula
        delta = J' \ (-F);
        v = v + delta(1:numel(v));
        lambda = lambda + delta(end);

        % Check for convergence
        if norm(delta) < tolerance
            disp(['Converged in ', num2str(iteration), ' iterations.']);
            break;
        end
    end

    % Display a warning if maximum iterations reached
    if iteration == max_iterations
        disp('Warning: Maximum iterations reached without convergence.');
    end
  
end


function F = Hosein_Nikkhah_HW2_P2_1_1(v, lambda, A)
    % Eigenvalue equation: Av - lambda*v
    eq1 = A * v - lambda * v;
    % Normalization constraint: ||v||^2 - 1
    eq2 = norm(v)^2 - 1;
    % Combine the equations into a vector
    F = [eq1; eq2];
end

function J = jacobian_matrix(v, lambda, A)
    % Assuming v is a column vector
    m = size(A, 1); % Number of rows in A
    n = numel(v) + 1; % Number of input variables (v + lambda)
    
    % Initialize Jacobian matrix
    J = zeros(m, n);
    
    % Compute partial derivatives
    for i = 1:m
        for j = 1:numel(v)
            J(i, j) = A(i, j) + lambda * (j == i);
            % Partial derivative of (Av)_i with respect to v_j
        end
        J(i, numel(v) + 1) = v(i);
        % Partial derivative of (lambda * v)_i with respect to lambda
    end
end


