H = 74;             % Height of jump point (m)
D = 31;             % Deck height (m)
c = 0.9;            % Drag coefficient (kg/m)
m = 80;             % Mass of the jumper (kg)
L = 43.6;             % Length of bungee cord (m)
k = 75.8;             % Spring constant of bungee cord (N/m)
g = 9.8;            % Gravitational acceleration (m/s^2)
C = c/m;            % Scaled drag coefficient
K = k/m;            % Scaled spring constant

T = 60;             % Final time in simulation (s)
n = 10000;          % Number of subintervals (equivalent to 0.006s / step)


[t, y, v, h] = rk4_bungee(T, n, g, C, K, L);

figure
plot(t, y);
xlabel('time (s)');
ylabel('distance fallen (m)');
title('Figure 1: Distance vs. Time');
max(y)
