Documentation of the solution from the Test_primal_active_set.m script

* Convex QP problem definition
H = [2 -1; -1 2];  % Hessian matrix (positive definite)
g = [-2; -6];      % Gradient vector
Aineq = [-1 1; 1 1];  % Inequality constraint matrix
bineq = [2; 10];     % Inequality constraint vector
lb = [0; 0];        % Lower bounds
ub = [inf; inf];    % Upper bounds


>> Test_primal_active_set

--- Solution Comparison ---

Primal Variables (x):
 quadprog:  x = [3.333333, 4.666667]
 Primal Active Set: x = [3.333333, 4.666667]

Dual Variables (Lagrange Multipliers):
  quadprog:  lambda = [0.000000, 0.000000]
  Active Set: z = [0.000000, 0.000000]

Objective Function Value:
  quadprog Objective value = -17.333333
  Active Set Objective Value = -17.333333