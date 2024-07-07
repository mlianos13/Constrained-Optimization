% Define QP Problem Data
H = [2 -1; -1 2];  % Hessian matrix (positive definite)
g = [-2; -6];      % Gradient vector
Aineq = [-1 1; 1 1];  % Inequality constraint matrix
bineq = [2; 10];     % Inequality constraint vector
lb = [0; 0];        % Lower bounds
ub = [5; 5];    % Upper bounds

%-----------------------------------------------------------------------------------------------------------------------------------

% Solve using quadprog (Standard QP Solver)
options_quadprog = optimoptions('quadprog', 'Display', 'off');
[x_quadprog, fval_quadprog, exitflag_quadprog, output, lambda] = quadprog(H, g, Aineq, bineq, [], [], lb, ub, [], options_quadprog);

%-----------------------------------------------------------------------------------------------------------------------------------

% Solve using the primal active-set algorithm 

% Find initial feasible point
n = size(H, 2);  % Number of original variables
x0 = findFeasibleInitialPoint(Aineq, bineq, lb, ub, n);

% Call the primal active set solver
tic; % Start timer
[x_opt_active_set, lambda_opt_active_set, active_set, iterations] = qpsolverActiveSet(H, g, Aineq', bineq, x0);
elapsed_time = toc; % Stop timer

% Calculate optimal objective value
optimal_obj_value = 0.5 * x_opt_active_set' * H * x_opt_active_set + g' * x_opt_active_set;

%---------------------------------------------------------------------------------------------------------------------------------

% Compare Solutions
fprintf('\n--- Solution Comparison ---\n');

fprintf('\nPrimal Variables (x):\n');
fprintf(' quadprog:  x = [%f, %f]\n', x_quadprog);
fprintf(' Primal Active Set: x = [%f, %f]\n', x_opt_active_set);

fprintf('\nDual Variables (Lagrange Multipliers):\n');
fprintf(' quadprog:  lambda = [%f, %f]\n', lambda.ineqlin);
fprintf(' Active Set: z = [%f, %f]\n', lambda_opt_active_set);

fprintf('\nNumber of Iterations:\n');
fprintf(' quadprog: n = %')
fprintf(' Active Set: n = %')

fprintf('\nObjective Function Value:\n');
fprintf(' quadprog Objective value = %f\n', fval_quadprog);
fprintf(' Active Set Objective Value = %f\n', optimal_obj_value);

