% Define your functions
syms x1 x2;
f1 = (x1 - 2)^2 + (x2 - 2)^2 - 9;
f2 = (x1 + 3)^2 + (x2 + 1)^2 - 9;

% Calculate the gradients (partial derivatives)
df1_x1 = diff(f1, x1);
df1_x2 = diff(f1, x2);
df2_x1 = diff(f2, x1);
df2_x2 = diff(f2, x2);

% Set up a grid for the contour plot
[x1_vals, x2_vals] = meshgrid(-10:0.1:10, -10:0.1:10);

% Calculate the values of your functions on the grid
f1_vals = double(subs(f1, {x1, x2}, {x1_vals, x2_vals}));
f2_vals = double(subs(f2, {x1, x2}, {x1_vals, x2_vals}));

% Create a contour plot
figure;
contour(x1_vals, x2_vals, f1_vals, [0, 0], 'DisplayName', 'f1(x1, x2)');
hold on;
contour(x1_vals, x2_vals, f2_vals, [0, 0], 'DisplayName', 'f2(x1, x2)');

% Initial guesses
x1_0 = 1;  % Initial guess for x1
x2_0 = 1;  % Initial guess for x2

% Stopping criteria and iteration limits
max_iterations = 100;
tolerance = 1e-6;

% Newton-Raphson iteration for the first function (f1)
x1 = x1_0;
x2 = x2_0;
for i = 1:max_iterations
    f1_x1_val = double(subs(df1_x1, {x1, x2}));
    f1_x2_val = double(subs(df1_x2, {x1, x2}));
    
    df1_x1_val = double(subs(f1, {x1, x2}));
    df1_x2_val = double(subs(f1, {x1, x2}));
    
    x1_new = x1 - f1_x1_val / df1_x1_val;
    x2_new = x2 - f1_x2_val / df1_x2_val;
    
    if norm([x1_new - x1, x2_new - x2]) < tolerance
        break;
    end
    
    % Plot the iteration
    plot([x1, x1_new], [x2, x2_new], 'r');
    
    x1 = x1_new;
    x2 = x2_new;
end

% Display the root for f1
fprintf('Root for f1: x1 = %f, x2 = %f\n', x1, x2);

% Repeat the same process for f2 if needed

% Add labels and legend
xlabel('x1');
ylabel('x2');
legend;

% Adjust the plot limits as needed
axis([-10 10 -10 10]);
