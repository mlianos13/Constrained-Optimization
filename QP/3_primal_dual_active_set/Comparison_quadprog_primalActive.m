load('C:\Users\mario\OneDrive\Υπολογιστής\MSc\Courses\2_Spring_Semester\1.)_Constrained_Optimization\Assignment 2024\matlab\QP_Test.mat');

% Define the inequality constraints (C'x <= du and -C'x <= -dl)
Aineq = [C'; -C'];
bineq = [du; -dl];

% Define the lower and upper bounds (l <= x <= u)
lb = l; 
ub = u;


%-----------------quadprog-------------------%

% Define the options for quadprog
options = optimoptions('quadprog',...
    'Display', 'iter-detailed', ...
    'Diagnostics','on', ...
    'Algorithm','interior-point-convex', ...
    'MaxIterations',200);

% Record the start time
tic;

% Calling quadprog to solve the optimization problem
[U, fval, exitflag, output, lambda] = quadprog(H, g, Aineq, bineq, [], [], [], [], [], options);

% Record the CPU time
cpuTime = toc;

% Display the results
fprintf('Number of iterations: %d\n', output.iterations);
fprintf('CPU time: %f seconds\n', cpuTime);

% More statistics
disp('More solution statistics:');
disp(output);

%----------------Primal Dual Active-Set-------------%

n = size(H);

% Combine all constraints into one matrix C_aug and one vector d_aug
C_aug = [eye(n); -eye(n); C'; -C'];
d_aug = [lb; -ub; dl; -du];

eta = 0.995;     % Step size parameter

% Compute initial point
n = size(H, 2);
[x0] = findInitialPointHeuristic(H, g, C_aug, d_aug, lb, ub);

% Solve the QP problem using the primal-dual interior-point method
tic; % Start timer
[x_opt, z_opt, s_opt, rL_, rC_, rSZ_] = primalDualFramework(H, g, C_aug', d_aug, x0, eta);
elapsed_time = toc; % Stop timer

% Calculate optimal objective value
optimal_obj_value = 0.5 * x_opt' * H * x_opt + g' * x_opt;

%----------------Compare Solutions--------------%


fprintf('\n--- Solution Comparison ---\n');

fprintf('\nPrimal Variables (x):\n');
fprintf(' quadprog:  x = [%f, %f]\n', U);
fprintf(' Primal Active Set: x = [%f, %f]\n', x_opt);

fprintf('\nDual Variables (Lagrange Multipliers):\n');
fprintf(' quadprog:  lambda = [%f, %f]\n', lambda.ineqlin);
fprintf(' Active Set: z = [%f, %f]\n', z_opt);

fprintf('\nNumber of Iterations:\n');
fprintf(' quadprog: n = %f\n', output.iterations)
fprintf(' Active Set: n = %f\n', 30)

fprintf('\nCPU Time:\n');
fprintf(' quadprog: n = %f\n', cpuTime)
fprintf(' Active Set: n = %f\n', elapsed_time)

fprintf('\nObjective Function Value:\n');
fprintf(' quadprog Objective value = %f\n', fval);
fprintf(' Active Set Objective Value = %f\n', optimal_obj_value);


norm_difference = norm(x_opt - U)
fprintf(' Norm of the Difference = %f\n', norm_difference);