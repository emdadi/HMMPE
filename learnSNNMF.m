function [Ohat, That, Pihat] = learnSNNMF(samples, X, Y)
    % ---------------------------------
    % samples = matrix where columns are sample-pairs from HMM
    % X       = number of hidden states
    % Y       = number of possible observations (Y >= X)
    %
    % OUTPUT:
    % Estimates of transition and emission matrix, as well as 
    % stationary distribution of HMM, where the columns sum to one
    % ---------------------------------
    
    % First estimate S_{2,1}:
    
    S_21 = zeros(Y, Y);
    
    for ii = 1 : length(samples)
        x1 = samples(1, ii);
        x2 = samples(2, ii);
        
        S_21(x2, x1) = S_21(x2, x1) + 1;
    end
    
    S_21 = S_21 / length(samples);

    % Formulate the constraints for the optimization problem:
    % x will be = [ Ohat(:); { T diag(pi) }(:) ]
    
    % Every element inside [0, 1]
    lb = zeros(Y * X + X * X, 1);
    ub =  ones(Y * X + X * X, 1);
    
    Aeq = [];
    
    % The sum of all elements of { T diag(pi) } should be one
    Aeq = [zeros(1, X*Y) ones(1, X*X)];
    
    % The columns of Ohat should sum to one
    for ii = 1 : X
        Aeq_row = [zeros(1, (ii - 1) * Y) ones(1, Y) zeros(1, (X - ii) * Y + X * X)];
        Aeq = [Aeq; Aeq_row];
    end
    
    beq = ones(1 + X, 1);
    
    % Now solve the optimization problem:
    
    % Initial guess:
    That = eye(X);
    Ohat = [eye(X); zeros(Y - X, X)];
    Pihat = ones(X, 1) / X;
    x0 = [Ohat(:); That(:) / X];
    
    % Use multistart (Z times) to try to find the global minima
    ms = MultiStart;
    problem = createOptimProblem('fmincon','x0',x0, ...
        'objective', @(x) value(S_21, X, Y, x), 'lb', lb, ...
        'ub', ub, 'Aeq', Aeq, 'beq', beq);
    
    ms_times = 25;
    x = run(ms, problem, ms_times);

    % Chop x to matrices
    A = zeros(X, X);
    for ii = 1 : X
        Ohat(:, ii) = x(1 + (ii - 1)*Y : ii * Y);
        A(:, ii) = x(1 + X*Y + (ii - 1) * X : X*Y + ii*X);
    end
    
    Pihat = A' * ones(X, 1);
    That  = A * inv(diag(Pihat));
end
