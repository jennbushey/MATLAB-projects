% ENGG 681
% Assignment 4 Problem 2
% Numerical Differentiation
% Instructor:  Sameh Nassar 
% Submission Date: April 3, 2024

t = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30];
L = [0, 102, 408, 948, 1538, 2376, 3100, 3983, 5744, 6500, 8409, 9407, 10700, 12200, 15095, 16121];
velocity = zeros(size(t));
acceleration = zeros(size(t));

% Find h
H = diff(t);
if sum(H)/length(H) == H(1)
    h = H(1);
    fprintf("h is equidistant = %d\n", h)
end


% Velocity Finite Difference with Truncation Error of O(h2)
% Velocity is the first derivative of altitude vs time

% Finite Difference Equations
firstD_3p_forward = @(f, h, i) (-3*f(i) + 4*f(i+1) - f(i+2)) / (2*h);
firstD_2p_central = @(f, h, i) (f(i+1) - f(i-1)) / (2*h);
firstD_3p_backward = @(f, h, i) (f(i-2) - 4*f(i-1) + 3*f(i)) / (2*h);

for i = 1:length(t)
    if i == 1
        velocity(i) = firstD_3p_forward(L, h, i);
    elseif i == length(t)
        velocity(i) = firstD_3p_backward(L, h, i);
    else
        velocity(i) = firstD_2p_central(L, h, i);
    end
end


% Accelearation Finite Difference with Truncation Error of O(h2)
% Acceleration is the second derivative of altitude vs time

% Finite Difference Equations
secondD_4p_forward = @(f, h, i) (2*f(i) - 5*f(i+1) +4*f(i+2) - f(i+3)) / h^2;
secondD_3p_central = @(f, h, i) (f(i-1) - 2*f(i) + f(i+1)) / h^2;
secondD_4p_backward = @(f, h, i) (-f(i-3) + 4 * f(i-2) - 5*f(i-1) + 2*f(i)) / h^2;

for i = 1:length(t)
    if i == 1
        acceleration(i) = secondD_4p_forward(L, h, i);
    elseif i == length(t)
        acceleration(i) = secondD_4p_backward(L, h, i);
    else
        acceleration(i) = secondD_3p_central(L, h, i);
    end
end

% Print data in a table format
fprintf('\nTime\tAltitude\tVelocity\tAcceleration\n');
for i = 1:length(t)
    fprintf('%d\t%f\t%f\t%f\n', t(i), L(i), velocity(i), acceleration(i));
end


% Plot altitude vs time
subplot(3,1,1);
plot(t, L, 'r-o'); 
xlabel('time');
ylabel('altitude');
title('Altitude vs Time');

% Plot velocity vs time
subplot(3,1,2);
plot(t, velocity, 'b-x'); 
xlabel('time');
ylabel('velocity');
title('Velocity vs Time');

% Plot acceleration vs time
subplot(3,1,3);
plot(t, acceleration, 'g-*'); 
xlabel('time');
ylabel('acceleration');
title('Acceleration vs Time');

% Add an overall title
sgtitle('Group 67: Altitude, Velocity, and Acceleration vs Time');


