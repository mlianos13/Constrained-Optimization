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

%----------------Primal Active-Set-------------%

% Find initial feasible point
n = size(H, 2);  % Number of original variables
x0 = findFeasibleInitialPoint(Aineq, bineq, lb, ub, n);

% Call the primal active set solver
tic; % Start timer
[x_opt_active_set, lambda_opt_active_set, active_set, iterations] = qpsolverActiveSet(H, g, Aineq', bineq, x0);
elapsed_time = toc; % Stop timer

% Calculate optimal objective value
optimal_obj_value = 0.5 * x_opt_active_set' * H * x_opt_active_set + g' * x_opt_active_set;


%----------------Compare Solutions--------------%


fprintf('\n--- Solution Comparison ---\n');

fprintf('\nPrimal Variables (x):\n');
fprintf(' quadprog:  x = [%f, %f]\n', U);
fprintf(' Primal Active Set: x = [%f, %f]\n', x_opt_active_set);

fprintf('\nDual Variables (Lagrange Multipliers):\n');
fprintf(' quadprog:  lambda = [%f, %f]\n', lambda.ineqlin);
fprintf(' Active Set: z = [%f, %f]\n', lambda_opt_active_set);

fprintf('\nNumber of Iterations:\n');
fprintf(' quadprog: n = %f\n', output.iterations)
fprintf(' Active Set: n = %f\n', iterations)

fprintf('\nCPU Time:\n');
fprintf(' quadprog: n = %f\n', cpuTime)
fprintf(' Active Set: n = %f\n', elapsed_time)

fprintf('\nObjective Function Value:\n');
fprintf(' quadprog Objective value = %f\n', fval);
fprintf(' Active Set Objective Value = %f\n', optimal_obj_value);



U_trimmed = U(10:end);
x_opt_active_set_trimmed = x_opt_active_set(10:end);



norm_difference = norm(x_opt_active_set_trimmed - U_trimmed)
fprintf(' Norm of the Difference = %f\n', norm_difference);




























% Plot settings
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15);

%% Load variables
load('plot_output.mat')

%% Simulate system and convert to physical variables
u = reshape(U,[2,100]);

x = x_kk;
z = C*x_kk;
for i = 1:N-1
    x(:,i+1) = A*x(:,i) + B*u(:,i);
    z(:,i+1) = C*x(:,i+1);
end

% Convert to physical variables
%x = x + [xs; ds];
z = z + zs;
u = u + us;
r = reshape(R_k,[2 100]) + zs;

%% Plot inputs
figure
tiledlayout(2,1)
nexttile
for i=1:2
    stairs(T,u(i,:),'Linewidth',2)
    hold on
end
title('Inputs')
xlabel('Time [min]')
ylabel('Flow rate [cm$^3$/s]')
yline(400)
yline(0)
ylim([-10,410])
xlim([T(1) T(end)])
legend('$F_1$','$F_2$','$u_{min}$','$u_{max}$','Location','southwest')

%% Plot system outputs

nexttile
hold on
for i = 1:2
    stairs(T,z(i,:),'Linewidth',2)
    stairs(T,r(i,:),'Linewidth',2,'LineStyle','--')
end
xlabel('Time [min]')
ylabel('Height [cm]')
title('Outputs')
xlim([T(1) T(end)])
legend('$z_1$','$r_1$','$z_2$','$r_2$')














% Plot settings
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15);

%% Load variables
load('plot_output.mat')

%% Simulate system and convert to physical variables
u = reshape(x_opt_active_set,[2,100]);

x = x_kk;
z = C*x_kk;
for i = 1:N-1
    x(:,i+1) = A*x(:,i) + B*u(:,i);
    z(:,i+1) = C*x(:,i+1);
end

% Convert to physical variables
%x = x + [xs; ds];
z = z + zs;
u = u + us;
r = reshape(R_k,[2 100]) + zs;

%% Plot inputs
figure
tiledlayout(2,1)
nexttile
for i=1:2
    stairs(T,u(i,:),'Linewidth',2)
    hold on
end
title('Inputs')
xlabel('Time [min]')
ylabel('Flow rate [cm$^3$/s]')
yline(400)
yline(0)
ylim([-10,410])
xlim([T(1) T(end)])
legend('$F_1$','$F_2$','$u_{min}$','$u_{max}$','Location','southwest')

%% Plot system outputs

nexttile
hold on
for i = 1:2
    stairs(T,z(i,:),'Linewidth',2)
    stairs(T,r(i,:),'Linewidth',2,'LineStyle','--')
end
xlabel('Time [min]')
ylabel('Height [cm]')
title('Outputs')
xlim([T(1) T(end)])
legend('$z_1$','$r_1$','$z_2$','$r_2$')