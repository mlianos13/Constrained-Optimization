% Define QP Problem Data
H = [2 -1; -1 2];  % Hessian matrix (positive definite)
g = [-2; -6];      % Gradient vector
Aineq = [-1 1; 1 1];  % Inequality constraint matrix
bineq = [2; 10];     % Inequality constraint vector
lb = [0; 0];        % Lower bounds
ub = [5; 5];        % Upper bounds

% Problem dimensions
n = size(H, 1);

% Combine all constraints into one matrix C_aug and one vector d_aug
C_aug = [eye(n); -eye(n); Aineq];
d_aug = [lb; -ub; bineq];

% Set parameters
eta = 0.995;     % Step size parameter

% Compute initial point
[x0] = findInitialPointHeuristic(H, g, C_aug, d_aug, lb, ub);

% Solve the QP problem using the primal-dual interior-point method
[x_opt, z_opt, s_opt, residuals] = primalDualFramework(H, g, C_aug', d_aug, x0, eta);


% Set parameters
maxIter = 30;   % Maximum number of iterations
tol = 1e-6;      % Tolerance for convergence
eta = 0.995;     % Step size parameter

% Compute initial point
[x0, z0, s0] = findInitialPointHeuristic(H, g, C_aug, d_aug, lb, ub);

% Solve the QP problem using the primal-dual interior-point method
[x_opt, z_opt, s_opt] = primalDualFramework(H, g, C_aug', d_aug, x0, eta);

% Display the optimal solution
disp('Optimal x:');
disp(x_opt);
disp('Optimal z:');
disp(z_opt);
disp('Optimal s:');
disp(s_opt);