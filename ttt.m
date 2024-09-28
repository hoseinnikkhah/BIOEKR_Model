% Read CSV file
data = readmatrix('data.csv');

API_csv = data(1, 1);       % Get the value for API
Mu_csv = data(1, 3);        % Get the value for viscosity
n_csv = data(1, 4);         % Get the value for Porosity
epsilon_csv = data(1, 6);   % Get the value for Relative Permittivity
tau_csv = data(1,9);        % Get the value for Tortuosity
growth_csv = data(1,10);    % Get the value for Microbial Growth
