% ENGG 681
% Assignment 3 Problem 2
% Interpolation
% Instructor:  Sameh Nassar 
% Submission Date: March 13, 2024


% ------------ DATA COMMON TO BOTH PROGRAMS ------------
% Define the data points
t=[20, 29, 45, 58, 77, 84, 105, 112];
u=[214, 275, 84, 129, 95, 116, 178, 103];

% Define the points to interpolate
t_interpolated = [50, 90];
u_interpolated=zeros(size(t_interpolated));



%
%
%
% ------------ LAGRANGE ------------
fprintf("\n Lagrange Interpolating Polynomials \n")
% Display the interpolated values
fprintf(' t     u_interpolated\n');
for i = 1:length(t_interpolated)
    % Call the lagrange_interpolation function (function at end of script)
    u_interpolated(i) = lagrange_interpolation(t, u, t_interpolated(i));
    % Print the interpolated values
    fprintf("%5.1f  %10.6f\n", t_interpolated(i), u_interpolated(i))
end



%
%
%
% ------------ QUADRATIC SPLINE ------------
fprintf("\n Quadratic Spline Interpolation \n")

n = length(t)-1;  % number of data points
k = 3*n;          % number of equations/unknowns

A = zeros(k);     % matrix of unknowns
b = zeros(k, 1);  % equation solutions

% The function values of adjacent polynomials must be equal at the interior knots.
for i = 2:n
    A(2 * (i - 1) - 1, 3 * (i - 1) - 2) = t(i)^2;
    A(2 * (i - 1) - 1, 3 * (i - 1) - 1) = t(i);
    A(2 * (i - 1) - 1, 3 * (i - 1) - 0) = 1;
    b(2 * (i - 1) - 1) = u(i);

    A(2 * (i - 1), 3 * (i - 1) + 1) = t(i)^2;
    A(2 * (i - 1), 3 * (i - 1) + 2) = t(i);
    A(2 * (i - 1), 3 * (i - 1) + 3) = 1;
    b(2 * (i - 1)) = u(i);
end

% The first and last functions must pass through the end points.
A(2 * n - 1, 1) = t(1)^2;       % a0
A(2 * n - 1, 2) = t(1);         % b0
A(2 * n - 1, 3) = 1;            % c0
b(2 * n - 1) = u(1);

A(2 * n , 3 * (n - 1) + 1) = t(n + 1)^2;    % an
A(2 * n , 3 * (n - 1) + 2) = t(n + 1);      % bn
A(2 * n , 3 * (n - 1) + 3) = 1;             % cn
b(2 * n ) = u(n + 1);

% The 1st derivatives at the interior knots must be equal
% Add conditions for smoothness at the interior knots
for i = 2:n
    A(2 * n + i - 1, 3 * (i - 1) - 2) = 2 * t(i);     % ai-1
    A(2 * n + i - 1, 3 * (i - 1) - 1) = 1;            % bi-1

    A(2 * n + i - 1, 3 * (i - 1) + 1) = -2 * t(i);    % ai-1
    A(2 * n + i - 1, 3 * (i - 1) + 2) = -1;           % bi-1
end


% Assume that the 2nd derivative is zero at the first point
A(3 * n, 1) = 2;

% View the matrices
% A
% b

coeffs = A \ b;

% View the coefficients
% coeffs

% Interpolate the values in t_interpolated
for i = 1:length(t_interpolated)
    for j = 2:length(t)
        if t_interpolated(i) <= t(j)
            u_interpolated(i) = ...
            coeffs(3 * (j-1) - 2) * t_interpolated(i)^2 + ...
            coeffs(3 * (j-1) - 1) * t_interpolated(i) + ...
            coeffs(3 * (j-1) - 0);          
            break
        end
    end
end


% Display the interpolated values
fprintf(' t     u_interpolated\n');
for i = 1:length(t_interpolated)
    % Print the interpolated values
    fprintf("%5.1f  %10.6f\n", t_interpolated(i), u_interpolated(i))
end


%
%
%
% ------------ FUNCTIONS ------------

function y_interp = lagrange_interpolation(x_data, y_data, x_interp)
    % Lagrange Interpolation
    n = length(x_data); % n is the length of the input array
    
    for i = 1:length(x_interp)
        x = x_interp(i);
        sum = 0;
        for j = 1:n
            L = 1;
            for k = 1:n
                if k ~= j
                    L = L * (x - x_data(k)) / (x_data(j) - x_data(k));
                end
            end
            sum = sum + y_data(j) * L;
        end
        y_interp = sum;
    end
end