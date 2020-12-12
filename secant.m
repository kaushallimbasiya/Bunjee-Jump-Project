function [x,i] = secant(f, x0, x1, tol, maxiters)
%SECANT Secant method
% [x,i] = secant(f, x0, x1, tol, maxiters) performs the secant
% method with f(x), starting at x_0 = x0 and x_1 = x1, and continuing
% until either |x_i+1 - x_i| <= tol, or maxiters iterations have
% been taken. The number of iterations, i, is also returned.
% An error is raised if the first input is not a function handle.
% A warning is raised if the maximum number of iterations is reached
% without achieving the tolerance.

if ~isa(f, 'function_handle')
    error('Your first input was not a function handle')
end

% Perform first iteration to fill x
i = 1;
x = x1 - f(x1)*((x1 - x0) / (f(x1) - f(x0)));

% Perform subsequent iterations, ending when tolerance met or iters exceed
% maximum
while abs(x - x1) > tol && i < maxiters
    % Increment counter
    i = i + 1;
    
    % Shift prev and current x values left
    x0 = x1;
    x1 = x;
    
    % Get next iteration value
    x = x1 - f(x1)*((x1 - x0) / (f(x1) - f(x0)));
end

% Check if difference still exceeds tolerance, this means loop hit maxiters
% before finding a suitable value
if abs(x - x1) > tol
    warning('Maximum number of iterations reached without achieving tolerance.')
end