% Define the parameters
L = 1; % Length of the domain
Nx = 100; % Number of grid points in x-direction
Ny = 3; % Number of components in the model
Nt = 1000; % Number of time steps
dx = L/(Nx-1); % Grid spacing in x-direction
dt = 0.001; % Time step size
T = ones(Nx,Ny); % Initial condition for T
C = zeros(Nx,Ny); % Initial condition for C
C(1,:) = [0.0449, 0.0283, 0.04]; % Boundary condition at x=0
sigma = sum([1,-1,1].*[1.5,0,-0.5].*[0.0449,0.0283,0.04]); % Sigma value
D = ones(Nx,Ny); % Diffusion coefficient
u = ones(Nx,Ny); % Velocity
v = ones(Nx,Ny); % Reaction rate
a = ones(Ny,Ny); % Coefficient matrix
lambda = ones(Nx,Ny); % Decay rate
Jw = 1; % Water flux
I = [0,0,0]; % Current density
z = [1,1,1]; % Charge number
F = 96485; % Faraday constant

% Define the boundary conditions
h0 = 0; % Transport boundary condition at x=0
hL = 0; % Transport boundary condition at x=L
E = zeros(Nx,1); % Boundary condition at x=L
E(end) = 0;

% Define the finite difference matrices
A = zeros(Nx,Nx);
B = zeros(Nx,Nx);
C_mat = zeros(Nx,Nx);
D_mat = zeros(Nx,Nx);
E_mat = zeros(Nx,Nx);

for i = 2:Nx-1
    A(i,i-1) = -D(i,1)*(1/dx^2);
    A(i,i) = (1/dt) + D(i,1)*(2/dx^2) + v(i,1);
    A(i,i+1) = -D(i,1)*(1/dx^2);
    
    B(i,i-1) = -D(i,2)*(1/dx^2);
    B(i,i) = (1/dt) + D(i,2)*(2/dx^2) + v(i,2);
    B(i,i+1) = -D(i,2)*(1/dx^2);
    
    C_mat(i,i-1) = -D(i,3)*(1/dx^2);
    C_mat(i,i) = (1/dt) + D(i,3)*(2/dx^2) + v(i,3);
    C_mat(i,i+1) = -D(i,3)*(1/dx^2);
    
    D_mat(i,i-1) = (1/dx)*F*z(1)*D(i,1);
    D_mat(i,i) = (1/dx)*F*z(2)*D(i,2);
    D_mat(i,i+1) = (1/dx)*F*z(3)*D(i,3);
    
    E_mat(i,i-1) = (1/dx)*F*z(1)*u(i,1)*C(i,1);
    E_mat(i,i) = (1/dx)*F*z(2)*u(i,2)*C(i,2);
    E_mat(i,i+1) = (1/dx)*F*z(3)*u(i,3)*C(i,3);
end

% Define the transport boundary conditions
A(1,1) = 1;
B(1,1) = 1;
C_mat(1,1) = 1;
D_mat(1,1) = F*z(1)*D(1,1)/dx;
D_mat(1,2) = F*z(2)*D(1,2)/dx;
D_mat(1,3) = F*z(3)*D(1,3)/dx;
E_mat(1,1) = F*z(1)*u(1,1)*C(1,1)/dx;
E_mat(1,2) = F*z(2)*u(1,2)*C(1,2)/dx;
E_mat(1,3) = F*z(3)*u(1,3)*C(1,3)/dx;

A(Nx,Nx) = 1;
B(Nx,Nx) = 1;
C_mat(Nx,Nx) = 1;
D_mat(Nx,Nx) = F*z(1)*D(Nx,1)/dx;
D_mat(Nx,Nx-1) = -F*z(2)*D(Nx,2)/dx;
D_mat(Nx,Nx-2) = F*z(3)*D(Nx,3)/dx;
E_mat(Nx,Nx) = F*z(1)*u(Nx,1)*C(Nx,1)/dx;
E_mat(Nx,Nx-1) = -F*z(2)*u(Nx,2)*C(Nx,2)/dx;
E_mat(Nx,Nx-2) = F*z(3)*u(Nx,3)*C(Nx,3)/dx;

% Define the time loop
for n = 1:Nt
    % Solve for C using the tridiagonal matrix algorithm
    b = zeros(Nx,1);
    b(1) = h0;
    b(Nx) = hL;
    for i = 2:Nx-1
        b(i) = T(i,1)*Jw + I(1)/(z(1)*F) + a(1,2)*I(2)/(z(2)*F) + a(1,3)*I(3)/(z(3)*F);
    end
    C(:,1) = tridiag(A,b);
    
    b = zeros(Nx,1);
    b(1) = h0;
    b(Nx) = hL;
    for i = 2:Nx-1
        b(i) = T(i,2)*Jw + I(2)/(z(2)*F) + a(2,1)*I(1)/(z(1)*F) + a(2,3)*I(3)/(z(3)*F);
    end
    C(:,2) = tridiag(B,b);
    
    b = zeros(Nx,1);
    b(1) = h0;
    b(Nx) = hL;
    for i = 2:Nx-1
        b(i) = T(i,3)*Jw + I(3)/(z(3)*F) + a(3,1)*I(1)/(z(1)*F) + a(3,2)*I(2)/(z(2)*F);
    end
    C(:,3) = tridiag(C_mat,b);
    
    % Solve for T using the tridiagonal matrix algorithm
    b = zeros(Nx,1);
    for i = 2:Nx-1
        b(i) = lambda(i,1)*T(i,1) + sigma*(E(i+1)-E(i))/dx;
    end
    T(:,1) = tridiag(D_mat,b);
    
    b = zeros(Nx,1);
    for i = 2:Nx-1
        b(i) = lambda(i,2)*T(i,2) + sigma*(E(i+1)-E(i))/dx;
    end
    T(:,2) = tridiag(D_mat,b);
    
    b = zeros(Nx,1);
    for i = 2:Nx-1
        b(i) = lambda(i,3)*T(i,3) + sigma*(E(i+1)-E(i))/dx;
    end
    T(:,3) = tridiag(D_mat,b);
    
    % Update the boundary condition at x=L
    E(end) = 0;
end

% Plot the results
x = linspace(0,L,Nx);
figure;
plot(x,T(:,1),'r',x,T(:,2),'g',x,T(:,3),'b');
xlabel('x');
ylabel('T');

% Define the tridiagonal matrix algorithm
function x = tridiag(A, b)
    n = length(b);
    x = zeros(n,1);
    for i = 2:n
        A(i,i-1) = A(i,i-1)/A(i-1,i-1);
        A(i,i) = A(i,i) - A(i,i-1)*A(i-1,i);
        b(i) = b(i) - A(i,i-1)*b(i-1);
    end
    x(n) = b(n)/A(n,n);
    for i = n-1:-1:1
        x(i) = (b(i) - A(i,i+1)*x(i+1))/A(i,i);
    end
end

