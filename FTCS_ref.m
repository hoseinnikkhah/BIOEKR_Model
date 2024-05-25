L = 0.02; % length of domain in x direction

I0 = 0.0000829*24*60*60; % Bq m^-2 day^-1
L1 = (0.693/53.2); % Decay constant day^-1
tmax = 120; % end time
nx = 90;% number of nodes in x direction
nt = 121; % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);
alpha= 2.5*10^-13*24*3600;
r = alpha*dt/dx^2; rr = 1 - 2*r-L1*dt;
v = 2.5*10^-6; % Convection velocity m day^-1
J0 = I0/sqrt(L1*alpha); % total inventory of Be-7 in soil
% --- Create arrays to save data for export
x = linspace(0,L,nx);
t = linspace(0,tmax,nt);
U = zeros(nx,nt);

% --- Set IC and BC
U(:,1)= 0;

% --- Loop over time steps

for m= 2:nt
    U(1,m) = J0; %--- Upper boundary
    U(end,m) = 0; %--- Lower boundary
    
    for i= 2:nx-1
        
        
        U(i,m) = r*U(i-1,m-1)+ rr*U(i,m-1)+ r*U(i+1,m-1);
        
        
    end
end

plot(t,U(10,:),'--');
xlabel('Time');
ylabel('depth(m)');
set(gca,'YDir','reverse','XAxisLocation','top');


figure(2)
fig = plot(U(:,1),x');
for k=2:tmax+1
    set(fig,'xdata',U(:,k),'ydata',x')
    set(gca,'YDir','reverse','XAxisLocation','top');
    title('Surface plot of solution.');
    ylabel('Distance (m)');
    pause(.1)
end
