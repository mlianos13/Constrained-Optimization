function [x, z, s, rL_, rC_, rSZ_] = primalDualFramework(H, g, C_aug, d_aug, x0, eta)
    % Ensure initial points are column vectors
    p = length(d_aug); 
    x = x0(:);
    s = ones(p, 1);
    z = ones(p, 1);
    max_iterations = 30;
    n = length(g);

    % Preallocate residuals for storage
    rL_ = zeros(n, max_iterations);
    rC_ = zeros(p, max_iterations);
    rSZ_ = zeros(p, max_iterations);

    S = diag(s);
    Z = diag(z);

    % Compute initial residuals
    rL = H * x + g - C_aug * z;
    rC = d_aug - C_aug' * x + s;
    rSZ = s .* z;
    mu = (z' * s) / p;

    % Parameters
    STOP = false;
    k = 0;
    tol = 1e-9; % Tolerance for convergence

    while ~STOP && k < max_iterations
        % Compute H_bar
        SZ = S \ Z;
        H_bar = H + C_aug * SZ * C_aug';

        % Symmetrize H_bar
        H_bar = (H_bar + H_bar') / 2;

        % Regularize H_bar
        reg_param = 1e-8;
        H_bar = H_bar + reg_param * eye(size(H_bar));

        % Check condition number
        cond_H_bar = cond(H_bar);
        fprintf('Iteration %d: Condition number of H_bar: %e\n', k, cond_H_bar);

        try
            % Perform Cholesky decomposition
            T = chol(H_bar, 'lower');
        catch
            warning('H_bar is not positive definite at iteration %d', k);
            return;
        end

        % Affine Direction
        rL_tilde = rL - C_aug * SZ * (rC - Z \ rSZ);
        rhs = -rL_tilde;

        delta_x_aff = T' \ (T \ rhs); 
        delta_z_aff = -SZ * (C_aug') * delta_x_aff + SZ * (rC - Z \ rSZ);
        delta_s_aff = -Z \ rSZ - (Z \ S) * delta_z_aff;

        % Compute the largest alpha_aff such that z + alpha_aff * delta_z_aff >= 0 and s + alpha_aff * delta_s_aff >= 0
        alphazs = (-[z; s] ./ [delta_z_aff; delta_s_aff]);
        alpha_aff = min([1; alphazs([delta_z_aff; delta_s_aff] < 0)]);

        % Compute the affine duality gap
        mu_aff = (z + alpha_aff * delta_z_aff)' * (s + alpha_aff * delta_s_aff) / p;

        % Compute the centering parameter
        sigma = (mu_aff / mu)^3;

        % Affine-Centering-Correction Direction
        rSZ_tilde = rSZ + delta_s_aff .* delta_z_aff - sigma * mu * ones(p, 1);
        rL_tilde = rL - C_aug * SZ * (rC - Z \ rSZ_tilde);
        
        rhs = -rL_tilde;
        
        delta_x_cor = T' \ (T \ rhs); 
        
        delta_z_cor = -SZ * (C_aug') * delta_x_cor + SZ * (rC - Z \ rSZ_tilde);
        delta_s_cor = -Z \ rSZ_tilde - (Z \ S) * delta_z_cor;

        % Step size
        alphazs = (-[z; s] ./ [delta_z_cor; delta_s_cor]);
        alpha = min([1; alphazs([delta_z_cor; delta_s_cor] < 0)]);

        % Ensure positivity of z and s
        alpha_ = eta * alpha;

        % Combine the directions
        x = x + alpha_ * delta_x_cor;
        z = max(z + alpha_ * delta_z_cor, 1e-10); % Ensure z remains positive
        s = max(s + alpha_ * delta_s_cor, 1e-10); % Ensure s remains positive

        % Update variables
        S = diag(s);
        Z = diag(z);

        % Recompute residuals
        rL = H * x + g - C_aug * z;
        rC = d_aug - C_aug' * x + s;
        rSZ = s .* z;
        mu = (s' * z) / p;

        % Store residuals
        rL_(:, k+1) = rL;
        rC_(:, k+1) = rC;
        rSZ_(:, k+1) = rSZ;

        k = k + 1;

        % Check convergence conditions
        fprintf('Iteration %d: norm(rL) = %e, norm(rC) = %e, mu = %e\n', k, norm(rL, inf), norm(rC, inf), mu);
        STOP = (norm(rL, inf) < tol) && (norm(rC, inf) < tol) && (abs(mu) < tol);
    end
    
    if k >= max_iterations
        warning('Maximum number of iterations reached without convergence.');
    end

    % Trim unused preallocated space
    rL_ = rL_(:, 1:k);
    rC_ = rC_(:, 1:k);
    rSZ_ = rSZ_(:, 1:k);
end
