L = 0.3; % length of domain in x direction [m]

I0 = (10^-8.2)*1000; 
R1 = (0.000000693/53.2); 
tmax = 120; % end time [days]

nx = 31;% number of nodes in x direction
nt = 172801; % number of time steps

dx = L/(nx-1);    % [m]
dt = tmax/(nt-1); % [day]
%dt is reported in day unit but it exactly equals to 1 minute
%might change it to seconds

alpha= (9.312*10^-9);           % [m^2/s]
% might need changing it to minutes



a = 69.76;                      % experimental parameter [mV]
b = -20.71;                     % experimental parameter
c = 0.15;                       % experimental parameter
T = 25+273;                     % Absolute temperature [K]
n = 0.52;                       % Porosity
tau = 0.80;                     % Tortuosity
z = 1;                          % Valency
rho = 1620*10^-6;               % Clay Density [kg/cm^3]
F = 96485.33212;                % Faraday constant [C/mol]
R = 8.314;                      % Gas constant [J/mol.K]
dEdx = 8/0.3;                   % Voltage gradient [V/m]
V = dEdx*L;                     % Voltage [V]
u_e = 0;                        % Electromigration mobility [m^2/V.day] 
zeta = -0.0027;                 % Zeta potential [V]
e_w = 1.7;                      % water permittivity 
mu_water = 1;                   % Water Viscosity [kg/(m.s)]
A = 0.0049;                     % Cross section [m^2]
Qf_r = 400;                     % flushing flow rate [ml/day]
Qf = (Qf_r*10^-6)/86400;        % flushing flow rate [m^3/s]
Va = 49 * 10^-9;                % Anode chamber volume [m^3]
Vc = Va;                        % Cathode chamber volume [m^3]

% Calculated constants

r = -alpha*tau*dt/dx^2; rr = 1 - 2*r-R1*dt;

D_star = n*tau*alpha;
u_star = (D_star*z*F)/(R*T);
k_eo = (-(e_w*zeta)/mu_water)*n*tau;

v = u_e + k_eo;                 % Convection velocity
beta = (-v*dEdx)/n;
beta_b = (v*dEdx)/Va;
% Dimensionless parameters

r1 = -D_star*dt/((dx^2)*n) ;rr1 = beta*dt/(2*dx); rrr = (R1*dt);
r1_b = D_star*dt/((dx^2)*Va); rr1_b = beta_b;

J0 = I0/sqrt(R1*alpha); % total inventory of Be-7 in soil
% --- Create arrays to save data for export

x = linspace(0,L,nx);
t = linspace(0,tmax,nt);
U = zeros(nx,nt);
S = zeros(nx,nt);

% --- Set IC and BC
U(:,1)= (10^-8.2)*1000;
S(:,1)= (10^-8.2)*1000;
S(:,2)= (10^-8.2)*1000;

% --- Loop over time steps

for m= 2:nt
    U(1,m) = J0; %--- Upper boundary
    U(end,m) = 0; %--- Lower boundary

    for i= 2:nx-1
        
        
        U(i,m) = r*U(i-1,m-1)+ rr*U(i,m-1)+ r*U(i+1,m-1);

        
    end
end


for m = 2:nt-1
    % Set boundary conditions for each time step

   
    %--- Upper boundary
    dcdx_left = (S(2, m) - S(1, m))/dx; 
    J_i_left = -D_star * dcdx_left - (V) * S(1, m) * dEdx;
    % Apply the boSndary condition to the left boundary for S
    S(1, m) = S(2, m) - (dt/Va)*(-J_i_left*A - Qf*S(1, m) + R1);

    %--- Lower boundary
    dcdx_right = (S(end,m) - S(nx-1))/dx;
    J_i_right = -D_star * dcdx_right - (v) * S(end, m) * dEdx;
    % Apply the boundary condition to the left boundary for S
    S(end, m) = S(nx-1, m) + (dt/Va)*(J_i_right*A - Qf*S(end, m) + R1);

    %S(1,m) = r1 * (S(2,m) - 2*S(1,m) + S(1,m-1)) + rr1 * (S(2,m) - S(1,m-1)) + rrr;
    for i = 2:nx-1

        S(i,m+1) = S(i,m) + r*(S(i+1,m) - 2*S(i,m) + S(i-1,m)) + rr1*(S(i+1,m) - S(i-1,m)) + rrr;
       %-1*U(i,m) = -U(i,m+1)  + r*(U(i+1,m) - 2*U(i,m) + U(i-1,m)) + rr1*(U(i+1,m) - U(i-1,m)) + rrr
        

    end
end

% --- Compare with exact solution at the end of the simulation

t1 = exp(-sqrt(R1/alpha)*x).*erfc((x./(2*sqrt(alpha*tmax)))-sqrt(R1*tmax));
t2 = exp(sqrt(R1/alpha)*x).*erfc((x./(2*sqrt(alpha*tmax)))+ sqrt(R1*tmax));
ue1 = (I0/(2*sqrt(R1*alpha)))*(t1-t2)+(I0/sqrt(R1*alpha))*exp(-sqrt(R1/alpha).*(x+0.002));

err= norm(U(:,nt)- ue1');

figure;
%plot(t,U(10,:),'--','DisplayName', 'U');

hold on;
plot(t,S(10,:),'--','DisplayName', 'S');

xlabel('Time (Days)');
ylabel('Conc (mg/kg)');

legend();