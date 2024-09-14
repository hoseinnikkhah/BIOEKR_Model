% Sample data points
P = [0, 0.01, 0.02, 0.07, 0.14, 0.23, 0.26];
B = [0,0.01,0.04,0.09,0.17,0.25,0.30];
R = [0,0.03,0.05,0.07,0.15,0.23,0.25];
Bl = [0,0.02,0.04,0.08,0.14,0.21,0.24];

% Sample data points
epsilon = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];

% Fit a second-degree polynomial (adjust the degree as needed)
p = polyfit(epsilon, P, 2);
b = polyfit(epsilon, B, 2);
r = polyfit(epsilon, R, 2);
bl = polyfit(epsilon, Bl, 2);

% Generate x-values for missing points or extrapolation (replace with desired range)
new_epsilon = linspace(min(epsilon), max(epsilon), 100);

% Use the fitted polynomial to estimate epsilon for new x-values
new_P = polyval(p, new_epsilon);
new_B = polyval(b, new_epsilon);
new_R = polyval(r, new_epsilon);
new_Bl = polyval(bl, new_epsilon);

% Plot original data and fitted curves
plot(epsilon, P, 'o', new_epsilon, new_P, 'm-', ...
     epsilon, B, 'o', new_epsilon, new_B, 'k-', ...
     epsilon, R, 'o', new_epsilon, new_R, 'r-', ...
     epsilon, Bl, 'o', new_epsilon, new_Bl, 'b-')

% Set legends
legend('Data Points (Model 1)', 'Model 1', ...
       'Data Points (Model 2)', 'Model 2', ...
       'Data Points (Model 3)', 'Model 3', ...
       'Data Points (Model 4)', 'Model 4')
