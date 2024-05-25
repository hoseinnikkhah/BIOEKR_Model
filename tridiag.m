function x = tridiag(A,b)
% Solves the tridiagonal system Ax = b using the Thomas algorithm
% Input: A - tridiagonal matrix (NxN)
%        b - right-hand side vector (Nx1)
% Output: x - solution vector (Nx1)

N = length(b);
x = zeros(N,1);

% Forward elimination
for i = 2:N
    m = A(i,i-1)/A(i-1,i-1);
    A(i,i) = A(i,i) - m*A(i-1,i);
    b(i) = b(i) - m*b(i-1);
end

% Backward substitution
x(N) = b(N)/A(N,N);
for i = N-1:-1:1
    x(i) = (b(i) - A(i,i+1)*x(i+1))/A(i,i);
end
end
