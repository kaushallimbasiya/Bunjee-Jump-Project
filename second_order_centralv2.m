function fdash = second_order_centeralv2(x0,x1,h)
%SECOND_ORDER_FORWARD Second order difference approximation for f'(x0).
%   FDASH = FIRST_ORDER_FORWARD(X0, X1, H) returns (X1 - X0) ./ (2*H);
%   Where X0 and X1 are numerically solved functions of f(x_n-1) and
%   f(x_n+1) when solving for f'(n).


% Use the second order centeral difference formula to approximate fdash.
fdash = ((x1) - (x0)) ./ (2*h);

end

