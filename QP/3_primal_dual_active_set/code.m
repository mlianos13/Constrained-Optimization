% Load data from .mat file
load('C:\Users\mario\OneDrive\Υπολογιστής\MSc\Courses\2_Spring_Semester\1.)_Constrained_Optimization\Assignment 2024\matlab\QP_Test.mat');

% Define the lower and upper bounds (l <= x <= u)
lb = l; 
ub = u; 
n = size(H);

% Combine all constraints into one matrix C_aug and one vector d_aug
C_aug = [eye(n); -eye(n); C'; -C'];
d_aug = [lb; -ub; dl; -du];

eta = 0.995;     % Step size parameter

% Compute initial point
n = size(H, 2);
[x0] = findInitialPointHeuristic(H, g, C_aug, d_aug, lb, ub);

% Solve the QP problem using the primal-dual interior-point method
[x_opt, z_opt, s_opt, rL_, rC_, rSZ_] = primalDualFramework(H, g, C_aug', d_aug, x0, eta);


%------------Plotting--------------%

% Plot settings
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15);

%% Load variables
load('plot_output.mat')

%% Simulate system and convert to physical variables
u = reshape(x_opt,[2,100]);

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

% Compute norms of residuals
norm_rL = vecnorm(rL_);
norm_rC = vecnorm(rC_);
norm_rSZ = vecnorm(rSZ_);

% Plot residuals
figure;
semilogy(1:length(norm_rL), norm_rL, '-o', 'DisplayName', 'norm(rL)');
hold on;
semilogy(1:length(norm_rC), norm_rC, '-s', 'DisplayName', 'norm(rC)');
semilogy(1:length(norm_rSZ), norm_rSZ, '-d', 'DisplayName', 'norm(rSZ)');
hold off;
xlabel('Iteration');
ylabel('Residual Norm');
legend;
title('KKT Residuals vs. Iterations');
grid on;