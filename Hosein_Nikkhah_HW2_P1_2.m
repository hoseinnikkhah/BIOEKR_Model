% Initial input vector
x = [1.5, 1.5, 2, 2]';
tolerance = 1e-6;
maxIterations = 100;
% Initialize arrays to store function values and x values for plotting
function_values = zeros(4, 5);
x_values = zeros(4, 5);

% Perform Newton's Method for the first 5 steps
for i = 1:5
    [x, F, ~, ~] = NewtonsMethod(x, tolerance, maxIterations);
    x_values(:, i) = x;
    function_values(:, i) = F;
end

% Create the first figure
figure;
for i = 1:4
    subplot(2, 2, i);
    plot(1:5, function_values(i, :), '-o');
    title(['f(X', num2str(i), ') vs. Iteration']);
    xlabel('Iteration');
    ylabel(['f(X', num2str(i), ')']);
end

% Create the second figure
figure;
for i = 1:4
    subplot(2, 2, i);
    plot(1:5, x_values(i, :), '-o');
    title(['X', num2str(i), ' vs. Iteration']);
    xlabel('Iteration');
    ylabel(['X', num2str(i)]);
end
