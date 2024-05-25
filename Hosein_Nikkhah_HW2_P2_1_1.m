function F = Hosein_Nikkhah_HW2_P2_1_1(v, lambda, A)
    % Eigenvalue equation: Av - lambda*v
    eq1 = A * v - lambda * v;
    % Normalization constraint: ||v||^2 - 1
    eq2 = norm(v)^2 - 1;
    % Combine the equations into a vector
    F = [eq1; eq2];
end