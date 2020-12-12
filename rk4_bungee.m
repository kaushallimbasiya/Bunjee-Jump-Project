function [t, y, v, h] = rk4_bungee(T, n, g, C, K, L)
% rk4_bungee Runge-Kutta Fourth Order method (RK4) for the bungee jumping
% model.
% [t, y, v, h] = rk4_bungee(T, n, g, C, K, L) performs RK4 on
% the bungee jumping model, taking n steps from t = 0 to t = T.
% The initial conditions are y(0) = 0 and v(0) = 0.
% The inputs g, C, K and L are parameters:
%   g = gravity (m/s^2)
%   C = drag coefficient (kg/m)
%   K = spring constant (N/m)
%   L = length of bungee rope (m)
% The outputs are the time array t, the solution arrays y and v, and the
% subinterval width h.

% Calculate step width
h = T/n;

% Generate time step iterations into a 1x(n+1) matrix
t = 0:h:T;

% Initialise solution arrays y and v
y = zeros(1,n+1);
v = zeros(1,n+1);

% Declare the two functions for use in RK4 k-iterations below
f1 = @(t,y,v) (g - C*abs(v)*v - max(0, K*(y - L)));
f2 = @(t,y,v) v;

% Iterate through all the time steps
for j = 1:n   
    % Calculate k1 for both v and y (slope at beginning of interval
    k1_v = h*f1(t(j),y(j),v(j));
    k1_y = h*f2(t(j),y(j),v(j));
    
    % Calculate k2 for both v and y (incremental slope on interval
    % midpoint)
    k2_v = h*f1(t(j) + h/2,y(j) + k1_y/2,v(j) + k1_v/2);
    k2_y = h*f2(t(j)+ h/2,y(j) + k1_y/2,v(j) + k1_v/2);
    
    % Calculate k3 for both v and y (further incremental slope on interval
    % midpoint)
    k3_v = h*f1(t(j) + h/2, y(j) + k2_y/2, v(j) + k2_v/2);
    k3_y = h*f2(t(j) + h/2, y(j) + k2_y/2, v(j) + k2_v/2);
    
    % Calculate k4 for both v and y (further incremental slope on interval
    % midpoint)
    k4_v = h*f1(t(j) + h, y(j) + k3_y, v(j) + k3_v);
    k4_y = h*f2(t(j) + h, y(j) + k3_y, v(j) + k3_v);
    
    % Perform the final RK4 calculation to get the current y and v iterates
    v(j+1) = v(j) + 1/6 * (k1_v+ 2*k2_v + 2*k3_v + k4_v);
    y(j+1) = y(j) + 1/6 * (k1_y+ 2*k2_y + 2*k3_y + k4_y);
end
end
