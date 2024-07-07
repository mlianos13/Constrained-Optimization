% Improved function to find initial feasible point
function x0 = findImprovedFeasibleInitialPoint(Aineq, bineq, lb, ub, n)
    % Ensure lb and ub are column vectors
    lb = lb(:);
    ub = ub(:);

    m = size(Aineq, 1);  % Number of constraints

    % Objective: Minimize the sum of slack variables
    f = [zeros(n, 1); ones(m, 1)];  % Slack variables have a cost of 1

    % Modify Aineq, b to include slack variables
    A_mod = [Aineq, -eye(m)];  % Add -I to handle slack variables for A*x <= b
    b_mod = bineq;  % Right hand side remains the same

    % Bounds for slack variables (non-negative)
    lb_slack = zeros(m, 1);  % Slack variables must be non-negative
    ub_slack = inf(m, 1);  % No upper bound on slack variables

    % Concatenate bounds for decision variables and slack variables
    lb_full = [lb; lb_slack];  % Lower bounds for decision variables and slack
    ub_full = [ub; ub_slack];  % Upper bounds for decision variables and slack

    % Options for linprog
    options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');  % Use interior-point method

    % Solve the linear program
    [x_full, ~, exitflag] = linprog(f, A_mod, b_mod, [], [], lb_full, ub_full, options);

    % Check if a feasible solution was found
    if exitflag ~= 1
        error('No feasible solution found.');
    end

    % Extract the original variables
    x0 = x_full(1:n);
end
