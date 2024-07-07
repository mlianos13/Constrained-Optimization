% Load data from .mat file
load('C:\Users\mario\OneDrive\Υπολογιστής\MSc\Courses\2_Spring_Semester\1.)_Constrained_Optimization\Assignment 2024\matlab\QP_Test.mat');


% Define the inequality constraints (C'x <= du and -C'x <= -dl)
Aineq = [C'; -C'];
bineq = [du; -dl];

% Define the lower and upper bounds (l <= x <= u)
lb = l; 
ub = u; 

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

% Extract and display the dual variables
fprintf('Dual variables (Lagrange multipliers):\n');
disp(lambda);



% Plotting


%--------------------------------------------------------------------------
%   Author:
%       Nicola Cantisani (nicca@dtu.dk)
%
%--------------------------------------------------------------------------
%   Description:
%       Reshapes and plots the open loop solution (inputs) of the optimal 
%       control problem for the 4-tank system (02612 assignment 2024).
%       The outputs of the system (level of the water in tank 1 and 2) 
%       are also plotted.
%
%   Inputs:
%       U       :     solution (column) vector (stacked inputs)
%
%--------------------------------------------------------------------------


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

