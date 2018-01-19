function out = value(S_21, X, Y, x)
    % Chop x to get O and A = T diag(pi)
    
    % x = [O(:); A(:)];
    
    O = zeros(Y, X);
    A = zeros(X, X);
    for ii = 1 : X
        O(:, ii) = x(1 + (ii - 1)*Y : ii * Y);
        A(:, ii) = x(1 + X*Y + (ii - 1) * X : X*Y + ii*X);
    end
    
    % Now optimize:
    % || S_21 - O A O' ||
    %
    % NOTE: Optimize sum of square of errors of elements
    % i.e. trace(B B') i.e. the Frobenius norm

    B = S_21 - O * A * O';
    out = ones(1, Y) * (B .* B) * ones(Y, 1);
end
