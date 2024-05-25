% Inputs
x1 = [1.5, 1.5, 2, 2];
x2 = [-1.5, 2.1, -1.9, 2.1, -1.9, 2.1, -1.5, 2.1];

tolerance_1 = 1e-4;
maxIterations_1 = 20;

tolerance_2 = 1e-6;
maxIterations_2 = 100;

% Solve for x1
[x1_solution, F1, J1, steps1] = NewtonsMethod(x1, tolerance_1, maxIterations_1);

% Solve for x2
[x2_solution, F2, J2, steps2] = NewtonsMethod(x2, tolerance_2, maxIterations_2);

% Display results
disp('Solution for x1:');
disp(x1_solution);
disp('Solution for x2:');
disp(x2_solution);

% Generate and display the table
disp('Table of steps for x1:');
disp('Iteration | Norm of x | Norm of Objective Value');
disp(steps1);

disp('Table of steps for x2:');
disp('Iteration | Norm of x | Norm of Objective Value');
disp(steps2);
