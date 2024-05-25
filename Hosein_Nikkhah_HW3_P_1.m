% Given matrices
H = [3, 1; 1, 2];
g = [-2; 1];

% Define the unconstrained optimization problem
fun = @(x) g.' * x + 0.5 * x.' * H * x;

% Solve the unconstrained problem
x_unconstrained = fminunc(fun, [0; 0]);

% Define the constrained optimization problem
constraint = @(x) x(1)^2 + x(2)^2 - 1;
nonlcon = @(x) deal([], constraint(x));

% Solve the constrained problem
options = optimoptions('fmincon', 'Algorithm', 'interior-point');
x_constrained = fmincon(fun, x_unconstrained, [], [], [], [], [], [], nonlcon, options);

% Display results
disp('Unconstrained minimum:')
disp(['x_unconstrained = ', num2str(x_unconstrained.')])
disp(['Cost at x_unconstrained: ', num2str(fun(x_unconstrained))])

disp('Constrained minimum:')
disp(['x_constrained = ', num2str(x_constrained.')])
disp(['Cost at x_constrained: ', num2str(fun(x_constrained))])
disp(['Constraint value at x_constrained: ', num2str(constraint(x_constrained))])
