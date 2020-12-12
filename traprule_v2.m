function I = traprule_v2(f,h) 
%TRAPRULE Trapezoidal rule integration. 
%   I = TRAPRULE(F, H) returns the trapezoidal rule approximation for 
%   the integral of f(x) from x=A to x=B, using N subintervals, 
%   Where 'f' is an array of numerically solved function values, and 'h' is
%   the interval size.


% Set first value to equal f(1).
S = f(1); 

% Iterate through variables from f(a+h) to f(b-h)
for j = 2:length(f)-1
    S = S + 2*f(j);
    
end

% Add on last term and multiply by h/2.
I = h/2 * (S + f(end));
end