function y = forward_eval(X, T, x)
%FORWARD_EVAL Evaluate Newton's forward difference form of the
%interpolating polynomial
% y = FORWARD_EVAL(X, T, x) returns y = Pn(x), where Pn is the
% interpolating polynomial constructed using the abscissas X and
% forward difference table T.

[m,n] = size(T);

if m ~= n
error('T must be square.');
end

% Set the initial y result to y0 (first column of T)
y = ones(size(x));
y = y .* T(1,1);

% Calculate polynomial components to build result in y
for k = 1:n-1
    
    % Calculate s for binomial coefficient
    h = X(k+1) - X(k);
    s = (x - X(1)) ./ h;
    
    % Initialise P to one because multiplying results
    P = ones(size(x));
    
    % Calculate numerator of binomial coefficient
    for i = 0:k-1
        P = P .* (s - i);
    end
    
    % Calculate final binomial coefficient
    P = P ./ factorial(k);
    
    % Multiply binomial coefficient by delta(k)y0
    P = P .* T(k+1,k+1);

    % Add the degree k component of the polynomial
    y = y + P;
end