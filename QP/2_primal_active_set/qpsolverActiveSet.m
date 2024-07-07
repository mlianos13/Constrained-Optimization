function [xopt, lambdaopt, Wset, it] = qpsolverActiveSet(H, g, A, b, x0)
    % QPSOLVERACTIVESET solves a convex QP with an active set algorithm
    %
    % Solves the qp
    %
    %   min 0.5 x' H x + g' x
    %    x
    %   s.t. A' x + b >= 0
    %
    % using a primal active set method and a feasible initial point, x0.
    %
    % Syntax: [xopt, lambdaopt, Wset, it] = qpsolverActiveSet(H, g, A, b, x0)

    % Tolerances hardcoded. Could be part of specifications
    tol = 1.0e-8;
    tolLx = tol;
    tolc = tol;
    tollambda = tol;
    tolp = tol;

    % Start with empty working set
    [n, m] = size(A);
    Wset = zeros(0, 1);
    IWset = (1:m)';
    lambda = zeros(m, 1);
    x = x0;

    % QP data
    gk = H * x + g;
    nablaxL = gk - A * lambda;
    c = A' * x + b;  % c(x) = A' x + b >= 0

    % Check if the initial point is optimal
    KKTstationarity = (norm(nablaxL, 'inf') < tolLx);
    KKTconditions = KKTstationarity; % the other conditions satisfied 

    % Main loop
    maxit = 100 * (n + m);
    it = 0;
    while (~KKTconditions && (it < maxit))
        it = it + 1;

        % Solve equality constrained QP
        Aw = A(:, Wset);
        bw = zeros(size(Aw, 2), 1);
        [p, lambdaWset] = EqualityQPsubproblem(H, gk, Aw, bw);

        if (norm(p, 'inf') > tolp) % p is non-zero
            % find binding constraint (if any)
            alpha = 1.0;
            idc = -1;
            nIWset = size(IWset, 1);
            for i = 1:nIWset
                pA = A(:, IWset(i))' * p;
                if pA < 0.0
                    alphapA = -c(IWset(i), 1) / pA;
                    if alphapA < alpha
                        alpha = alphapA;
                        idc = i;
                    end
                end
            end
            % Take step, update data and working set
            x = x + alpha * p;
            gk = H * x + g;
            c = A' * x + b;    
            if idc > 0
                Wset = [Wset; IWset(idc)];
                IWset = [IWset(1:idc-1); IWset(idc+1:end)];
            end
        else % p is zero
            % find minimum lambda
            idlambdaWset = -1;
            minlambdaWset = 0.0;
            nWset = size(Wset, 1);
            for i = 1:nWset
                if lambdaWset(i) < minlambdaWset 
                    idlambdaWset = i;
                    minlambdaWset = lambdaWset(i);
                end
            end

            if idlambdaWset > 0 % update the working set, x = x
                % if minimum lambda < 0 remove constraint from working set
                IWset = [IWset; Wset(idlambdaWset, 1)];
                Wset = [Wset(1:idlambdaWset-1, 1); Wset(idlambdaWset+1:end, 1)];
            else % optimal solution found
                KKTconditions = true;
                xopt = x;
                lambdaopt = zeros(m, 1);
                lambdaopt(Wset, 1) = lambdaWset;
            end
        end
    end

    if ~KKTconditions
        xopt = [];
        lambdaopt = [];
        Wset = [];
    end
end
