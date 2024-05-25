function [eigenvalue, eigenvector] = find_eigenpair(A)
    % Initialize variables
    v = randn(size(A, 1), 1);  % Initial random eigenvector guess
    lambda = randn();  % Initial random eigenvalue guess
    tolerance = 1e-8;
    max_iterations = 1000;

    % Iterative Newton-Raphson process
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
        if norm(F) < tolerance
            break;
        end
    end

    eigenvalue = lambda;
    eigenvector = v;
    
    % Check against MATLAB's eig function
    [matlab_eigenvectors, matlab_eigenvalues] = eig(A);
    [~, index] = min(abs(diag(matlab_eigenvalues) - eigenvalue));
    matlab_eigenvector = matlab_eigenvectors(:, index);

    % Display results
    disp('Newton-Raphson Method Result:');
    disp('Eigenvalue:');
    disp(eigenvalue);
    disp('Eigenvector:');
    disp(eigenvector);

    disp('MATLAB eig Function Result:');
    disp('Eigenvalue:');
    disp(matlab_eigenvalues(index, index));
    disp('Eigenvector:');
    disp(matlab_eigenvector);
end

function J = jacobian_matrix(v, lambda, A)
    m = size(A, 1);
    n = numel(v) + 1;
    J = zeros(m, n);
    
    for i = 1:m
        for j = 1:numel(v)
            J(i, j) = A(i, j) + lambda * (j == i);
        end
        J(i, numel(v) + 1) = v(i);
    end
end

