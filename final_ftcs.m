L = 40; % length of domain in x direction [cm]

I0 = 0.0000829*24*60*60; % Bq m^-2 day^-1
R1 = (0.693/53.2); % Decay constant day^-1
tmax = 35; % end time [days]

nx = 41;% number of nodes in x direction
nt = 50401; % number of time steps

dx = L/(nx-1);    %[cm]
dt = tmax/(nt-1); %[day]

alpha= (3.5447*10^-9)*24*3600*10^4;    %[cm^2/day]
r = alpha*dt/dx^2; rr = 1 - 2*r-R1*dt;

T = 25+273;                            % Absolute temperature [K]
mu_oil = 51*36*24;                     % Oil Viscosity [kg/(cm.day)]
mu_water = 0.001*36*24;                % Water Viscosity [kg/(cm.day)]

dEdx = 1.5;                       % Voltage gradient [V/cm]
V = dEdx*L;                       % Voltage [V]
u_e = 0;                          % Electromigration mobility [m^2/V.day] 
zeta = -0.0027;                   % Zeta potential [V]
e_0 = 8.854*10^-14;               % permittivity of free space [F/cm]
e_r = 7.5/100;                    % relative permittivity of clay [F/cm]
e_oil = 2.3/100;                  % Crude oil permittivity [F/cm]
e_clay = e_0*e_r;                 % Clay permittivity [F/cm]
e_w = 0.05;                       % Water permittivity [F/cm]
n = 0.64;                         % Porosity
tau = 0.44;                       % Tortuosity
z = 1;                            % Valency
rho = 1620*10^-6;                 % Clay Density [kg/cm^3]

% Calculated constants
k_eo1 = (e_oil*zeta)*n/(mu_oil);                 % Electroosmotic mobility Alshawabkeh [m^2/V.s]
k_eo2 = (e_oil*zeta)*n/(mu_oil*(tau^2));         % Electroosmotic mobility Vane [m^2/V.s]
k_eo3 = (e_oil*zeta*dEdx)/(mu_oil);              % Electroosmotic mobility Shapiro [S/m]

v0 = 2.5*10^-6; % Convection velocity m day^-1

v = u_e + k_eo3; % Convection velocity m day^-1

beta = v*dEdx;

rr1 = beta*dt/(2*dx); rrr = (R1*dt);


J0 = I0/sqrt(R1*alpha); % total inventory of Be-7 in soil

% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);
U = zeros(nx,nt);
S = zeros(nx,nt);
X = zeros(nx,nt);
R = zeros(nx,nt);
E = zeros(nx,nt);
W = zeros(nx,nt);
II = zeros(nx,nt);
% --- Set IC and BC
U(:,1)= 2000;
S(:,1)= 12;
S(:,2)= 12;
X(:,1)= 4000;
R(:,1)= 6000;
E(:,1)= 8000;
W(:,1)= 10000;
II(:,1)= 2000;
% --- Loop over time steps

for m= 2:nt
    U(1,m) = J0; %--- Upper boundary
    U(end,m) = 0; %--- Lower boundary
    X(1,m) = J0; %--- Upper boundary
    X(end,m) = 0; %--- Lower boundary
    R(1,m) = J0; %--- Upper boundary
    R(end,m) = 0; %--- Lower boundary
    E(1,m) = J0; %--- Upper boundary
    E(end,m) = 0; %--- Lower boundary    
    W(1,m) = J0; %--- Upper boundary
    W(end,m) = 0; %--- Lower boundary      
    for i= 2:nx-1
        
        
        U(i,m) = r*U(i-1,m-1)+ rr*U(i,m-1)+ r*U(i+1,m-1);
        X(i,m) = r*X(i-1,m-1)+ rr*X(i,m-1)+ r*X(i+1,m-1);
        R(i,m) = r*R(i-1,m-1)+ rr*R(i,m-1)+ r*R(i+1,m-1);
        E(i,m) = r*E(i-1,m-1)+ rr*E(i,m-1)+ r*E(i+1,m-1);
        W(i,m) = r*W(i-1,m-1)+ rr*W(i,m-1)+ r*W(i+1,m-1);
        
    end
end

for m= 2:nt-1
   
    II(1,m) = J0; %--- Upper boundary
    II(end,m) = 0; %--- Lower boundary    
    for i= 2:nx-1
        
        II(i,m) = II(i,m+1) - r*(II(i+1,m) - 2*II(i,m) + II(i-1,m)) - rr*(II(i+1,m) + II(i-1,m)) - rrr;
        
    end
end

for m = 2:nt-1
    % Set boundary conditions for each time step
    S(1,m) = J0; %--- Upper boundary
    S(end,m) = 0; %--- Lower boundary
    
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

%plot(t,X(10,:),'--','DisplayName', 'X');

%plot(t,R(10,:),'--','DisplayName', 'R');

%plot(t,E(10,:),'--','DisplayName', 'E');

%plot(t,W(10,:),'--','DisplayName', 'W');

xlabel('Time (Days)');
ylabel('Conc (mg/kg)');
set(gca,'YDir','reverse','XAxisLocation','top');
legend();

%figure(2)
%fig = plot(U(:,1),x');
%for k=2:tmax+1
%    set(fig,'xdata',U(:,k),'ydata',x')
%    set(gca,'YDir','reverse','XAxisLocation','top');
%    title('Surface plot of solution.');
%    ylabel('Distance (m)');
%    pause(.1)
%end

%figure(3)
%fig = plot(S(:,1),x');
%for k=2:tmax+1
%    set(fig,'xdata',S(:,k),'ydata',x')
%    set(gca,'YDir','reverse','XAxisLocation','top');
%    title('Surface plot of solution of new FTCS');
%    ylabel('Distance (m)');
%    pause(.1)
%end

