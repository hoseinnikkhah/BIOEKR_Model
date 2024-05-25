function [x, F, J, steps] = NewtonsMethod(x0, tol, maxIter)
    x = x0;
    steps = [];
    
    for iter = 1:maxIter
        [F, J] = FunctionEvaluator(x);
        A = det(J);
        if A == 0
             disp ('Jacobin is singular')
             break
        end        
        dx = -J \ F;  
        x = x + dx; 
        
        % Calculate the norms
        xNorm = norm(x);
        FNorm = norm(F);
        
        % Store the step information
        steps = [steps; [iter, xNorm, FNorm]];
        
        % Check for convergence
        if FNorm < tol
            break;  % if Converged
        end
    end
    
end


