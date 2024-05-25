function [v_new, lambda_new] = Hosein_Nikkhah_HW2_P2_1_2(v, lambda, A)
    % Calculate the Jacobian matrix
    J = jacobian_matrix(v, lambda, A);
    % Calculate the function values
    F = Hosein_Nikkhah_HW2_P2_1_1(v, lambda, A);
    % Update using the Newton-Raphson formula
    delta = J' \ (-F);
    v_new = v + delta(1:numel(v));
    lambda_new = lambda + delta(end);
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
            J(i, j) = A(i, j) + lambda * (j == i); % Partial derivative of (Av)_i with respect to v_j
        end
        J(i, numel(v) + 1) = v(i); % Partial derivative of (lambda * v)_i with respect to lambda
    end
end
