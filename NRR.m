function [root, steps] = NewtonRaphson(func, func_derivative, initial_guess, tol, max_iterations)
    % Initialize variables
    x = initial_guess;
    steps = [];
    
    for iter = 1:max_iterations
        f = func(x);
        f_prime = func_derivative(x);
        
        % Calculate the next estimate using the Newton-Raphson formula
        x_next = x - f / f_prime;
        
        % Calculate the absolute error to check for convergence
        error = abs(x_next - x);
        
        % Store iteration information
        steps = [steps; [iter, x, f, error]];
        
        % Check for convergence
        if error < tol
            root = x_next;
            return;
        end
        
        % Update the current estimate
        x = x_next;
    end
    
    % If max iterations reached without convergence, return the last estimate
    root = x;
end

% Define the function and its derivative
func = @(x) x^3 - 2*x - 5;
func_derivative = @(x) 3*x^2 - 2;

% Initial guess, tolerance, and maximum iterations
initial_guess = 2;  % Choose an appropriate initial guess
tolerance = 1e-6;
max_iterations = 100;

% Call the Newton-Raphson method
[root, steps] = NewtonRaphson(func, func_derivative, initial_guess, tolerance, max_iterations);

% Display the results
fprintf('Root found: %f\n', root);
fprintf('Number of iterations: %d\n', steps(end, 1));

% Display a table of steps
disp('Table of Steps:');
disp('Iteration | x_i     | f(x_i)  | Absolute Error');
disp(steps);
