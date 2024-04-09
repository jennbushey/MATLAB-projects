% ENGG 681
% Assignment 1 Problem 2
% Numerical Methods for Determining Roots
% Instructor:  Sameh Nassar 
% Submission Date: January 31, 2024

% Bisection Method
fprintf("\n***Bisection Method***\n")

% Define the function
P = @(m) (-m + 4) * exp(-0.5 * m) - 2;

% Initial guesses
ml = 0.5;
mu = 1.5;
error = 0.005;

% Check if the initial guesses have opposite signs
if P(ml) * P(mu) > 0
    fprintf('Initial guesses must have opposite signs.');
else
    % Initialize variables
    iteration = 0;
    mrprev = Inf;
    apr = Inf;   % Initialize apr with a value greater than 0.005
    approx_root = Inf; % initialize root with a value of infinity
    
    % Print table headings
    fprintf('Iteration     ml          mu          mr        P(ml)      P(mr)      P(ml)*P(mr)   approximate %% relative error\n');

    % Bisection method loop
    while apr >= error
        mr = (ml + mu) / 2;
        
        % approximate relative error
        apr = abs((mr - mrprev) / mr) * 100;
        
        % Print the results
        fprintf('%7i %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %9.3f%%\n', iteration, ml, mu, mr, P(ml), P(mr), P(ml)*P(mr) ,apr);
    
        % Update interval [ml, mu]
        if P(mr) * P(ml) < 0
            mu = mr;
        elseif P(mr) * P(ml) > 0 
            ml = mr;
        else
            break;
        end
        mrprev = mr;
        approx_root = mr;
        iteration = iteration + 1;
    end
end

 fprintf("Approximate root is: %.6f\n",approx_root);

%%
% Newton-Raphson
fprintf("\n***Newton-Raphson***\n")

% Define the function
P = @(m) (-m + 4) * exp(-0.5 * m) - 2;
Pprime = @(m) exp(-0.5 * m)*(0.5*m - 3); 

%P = @(x) 0.95*x^3 - 5.9*x^2 + 10.9*x -6;
%Pprime = @(x) 3*0.95*x^2 - 2*5.9*x + 10.9;

% Initial guess
m0 = 1.5;
error = 0.005;

% Initialize variables
iteration = 0;
apr = Inf;   % Initialize apr with a value greater than 0.005
approx_root = Inf; % initialize root with a value of infinity

% Print table headings
fprintf("Iteration     m0        P(m0)      P'(m0)     approximate %% relative error\n");

% Print iteration 0
fprintf('%7i %11.6f %11.6f %11.6f %11.3f%%\n', iteration, m0, P(m0), Pprime(m0), apr);

% Newton-Raphson method loop
while apr >= error
    
    % general iteration formula
    m1 = m0 - (P(m0) / Pprime(m0));
      
    % approximate relative error
    apr = abs((m1 - m0) / m1) * 100;

    % set variables for next iteration
    m0 = m1;
    iteration = iteration + 1;

    approx_root = m0;
    
    % Print the results
    fprintf('%7i %11.6f %11.6f %11.6f %11.3f%%\n', iteration, m0, P(m0), Pprime(m0), apr);
end

fprintf("Approximate root is: %.6f\n",approx_root);

%%
% Secant
fprintf("\n***Secant***\n")

% Initial guess
m0 = 0.7;
m_1 = 0.2; % m minus 1
error = 0.005;

% Define the function
P = @(m) (-m + 4) * exp(-0.5 * m) - 2;

% Initialize variables
iteration = 0;
apr = Inf;   % Initialize apr with a value greater than 0.005
approx_root = Inf; % initialize root with a value of infinity

% Print table headings
fprintf("Iteration     m-1      P(m-1)        m0         P(m0)     approximate %% relative error\n");

% Print iteration 0
    fprintf('%7i %11.6f %11.6f %11.6f  %11.6f %9.3f%%\n', iteration, m_1, P(m_1), m0, P(m0), apr);

% Secant method loop
while apr >= error
    % iteration formula
    m1 = m0 - P(m0) * ((m0 - m_1) / (P(m0) - P(m_1)));

    % approximate relative error
    apr = abs((m1 - m0) / m1) * 100;

    % set variables for next iteration
    m_1 = m0;
    m0 = m1;
    iteration = iteration + 1;
    approx_root = m0;
    % Print the results
    fprintf('%7i %11.6f %11.6f %11.6f  %11.6f %9.3f%%\n', iteration, m_1, P(m_1), m0, P(m0), apr);
end

 fprintf("Approximate root is: %.6f\n",approx_root);