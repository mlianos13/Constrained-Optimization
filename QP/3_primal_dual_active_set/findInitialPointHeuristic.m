function [x0, z0, s0] = findInitialPointHeuristic(H, g, C_aug, d_aug, lb, ub)
    % Number of decision variables
    n = size(H, 2);

    % Initialize x_bar as the midpoint of the bounds
    lb = lb(:);
    ub = ub(:);
    x_bar = (lb + ub) / 2;

    % Initialize z_bar and s_bar to be small positive constants
    z_bar = ones(length(d_aug), 1);
    s_bar = ones(length(d_aug), 1);

    % Ensure z_bar and s_bar are positive
    if any(z_bar <= 0) || any(s_bar <= 0)
        error('Initial values of z_bar and s_bar must be positive.');
    end

    % Compute residuals
    rL = H * x_bar + g - C_aug' * z_bar;
    rC = d_aug - C_aug * x_bar + s_bar;
    rSZ = s_bar .* z_bar;

    % Compute H_bar and LDL factorization
    S_inv_Z = diag(s_bar ./ z_bar);
    H_bar = H + C_aug' * S_inv_Z * C_aug;
    [L, D, P] = ldl(H_bar);

    % Compute affine search direction
    rL_tilde = rL - C_aug' * S_inv_Z * (rC - diag(1 ./ z_bar) * rSZ);
    delta_x_aff = P * (D \ (L \ (P' * rL_tilde)));  % Solve H_bar * delta_x_aff = rL_tilde
    delta_z_aff = S_inv_Z * (rC - C_aug * delta_x_aff);
    delta_s_aff = -diag(1 ./ z_bar) * rSZ - diag(1 ./ z_bar) * diag(s_bar ./ z_bar) * delta_z_aff;

    % Compute initial point
    x0 = x_bar + delta_x_aff;
    z0 = max(1, z_bar + delta_z_aff);
    s0 = max(1, s_bar + delta_s_aff);
end