function [F, J, sparsity_pattern, bandwidth] = Hosein_Nikkhah_HW2_P1_3(x)
    n = length(x);  
    if n == 0, error('Input vector, x, is empty.'); end
    if mod(n, 2) ~= 0 
       error('Input vector, x ,must have an even number of components.');
    end
    
    odds  = 1:2:n;
    evens = 2:2:n;
    F = zeros(n, 1);
    
    F(odds, 1)  = (1 - x(odds)).^2;
    F(evens, 1) = 10 .* ((x(odds) - x(evens).^2)).^2; 
    
    % Evaluate the Jacobian matrix if nargout > 1 
    c = 2 .* (x(odds) - 1);
    d = 40 .* (x(evens) .* (x(evens).^2 - x(odds)));
    e = 20 .* (x(odds) - x(evens).^2);
    
    row_indices = [odds, evens, evens];
    col_indices = [odds, evens, odds];
    values = [c; d; e];
    J = sparse(row_indices, col_indices, values, n, n);
        
    sparsity_pattern = spones(J);
    bandwidth = max(sum(sparsity_pattern)) - 1; 
    disp('Sparsity Pattern:');
    disp(sparsity_pattern);
    disp(['Bandwidth: ', num2str(bandwidth)]);

end
