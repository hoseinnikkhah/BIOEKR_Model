%ODE solver with initial points

%This is a genral solver for ODE problems

%This solver solves this problem(LaTex): x^{2}\frac{d^{3} y}{dx^{3}} -2x\frac{d y}{dx} + 5y = 0
%Go to https://latexeditor.lagrida.com/ to check the problem

%this problem is 3rd order ODE that can't be solved with that order since
%matlab only solves first order problems, so reducing orders is needed to
%do so \frac{d y}{dx} is equal to z and \frac{d^{3} y}{dx^{3}} is equal to
% \frac{d^{2} z}{dx^{2}} but it still need reducing so there's another
% value called Z, \frac{d^{2} z}{dx^{2}} equals to Z

function ODESolverScript()
    % Define the odesolver function
    format('compact')
    clc;
    function p = odesolver(x1, q)
        p(1, 1) = q(2); %q(2) is z in reduced form for \frac{d y}{dx}
        p(2, 1) = q(3); %q(2) is Z in reduced form for \frac{d^{3} y}{dx^{3}}
        p(3, 1) = 5 * q(1) / x1^2 - 2 * q(2) / x1; %q(1) is y in original form of problem
    end

    % Call ode23 within the script
    [x, u] = ode23(@odesolver, [1, 9], [2, 1, 0]);

    % Plot the results
    plot(x, u);
    xlabel('(x)');
    ylabel('(u)');
end
