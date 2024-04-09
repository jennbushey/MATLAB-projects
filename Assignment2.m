% ENGG 681
% Assignment 2 Problem 2
% Linear Equations and Systems
% Instructor:  Sameh Nassar 
% Submission Date: February 14, 2024

%%%%%% Define the functions
% Original Functions
% F1 = @(j, k, l, m, n) -1*j - 2*k + 4*l - 2*m + 0*n - 97.54; % = 0
% F2 = @(j, k, l, m, n)  1*j - 3*k + 2*l + 0*m - 1*n - 79.75; % = 0
% F3 = @(j, k, l, m, n)  0*j + 2*k - 3*l + 6*m + 0*n - 102.28; % = 0
% F4 = @(j, k, l, m, n)  3*j + 0*k + 1*l - 2*m + 7*n - 84.49; % = 0
% F5 = @(j, k, l, m, n) -3*j - 1*k + 0*l + 1*m - 2*n - 102.37; % = 0

% Re-arrange functions to be diagnonally dominant
F5 = @(j, k, l, m, n) -3*j - 1*k + 0*l + 1*m - 2*n - 102.37; % = 0
F2 = @(j, k, l, m, n)  1*j - 3*k + 2*l + 0*m - 1*n - 79.75; % = 0
F1 = @(j, k, l, m, n) -1*j - 2*k + 4*l - 2*m + 0*n - 97.54; % = 0
F3 = @(j, k, l, m, n)  0*j + 2*k - 3*l + 6*m + 0*n - 102.28; % = 0
F4 = @(j, k, l, m, n)  3*j + 0*k + 1*l - 2*m + 7*n - 84.49; % = 0

% Solve for the diagonally dominant variable
Fj = @(j, k, l, m, n) (102.37 - (-1*k + 0*l + 1*m - 2*n)) / -3;
Fk = @(j, k, l, m, n) (79.75  - ( 1*j + 2*l + 0*m - 1*n)) / -3;
Fl = @(j, k, l, m, n) (97.54  - (-1*j - 2*k - 2*m + 0*n)) / 4;
Fm = @(j, k, l, m, n) (102.28 - ( 0*j + 2*k - 3*l + 0*n)) / 6;
Fn = @(j, k, l, m, n) (84.49  - ( 3*j + 0*k + 1*l - 2*m)) / 7;

% Approximation Error Formula
% Ea = current approximation - previous approximation
Ferror = @(old, new) abs((new - old));



%%%%%% Gauss-Seidel
fprintf("\nGauss Seidel Method\n")

% Initial guesses
j = 0;
k = 0;
l = 0;
m = 0;
n = 0;
error = 0.00001; % approximation error

% Initialize Variables
iteration = 1;
approx_error = Inf;

% Begin Iteration
fprintf("iter  unknown  value       error       maximum error \n");
while approx_error > error
    % Calculate new values using most recent value
    j_new = Fj(j, k, l, m, n);
    k_new = Fk(j_new, k, l, m, n);
    l_new = Fl(j_new, k_new, l, m, n);
    m_new = Fm(j_new, k_new, l_new, m, n);
    n_new = Fn(j_new, k_new, l_new, m_new, n);
    
    % Calculate error for each variable
    error_j = Ferror(j, j_new);
    error_k = Ferror(k, k_new);
    error_l = Ferror(l, l_new);
    error_m = Ferror(m, m_new);
    error_n = Ferror(n, n_new);

    % Continue iteration based on maximum error the set of variables
    approx_error = max([error_j, error_k, error_l, error_m, error_n]);

    j = j_new;    
    k = k_new;
    l = l_new;
    m = m_new;
    n = n_new;

    fprintf("%3i      j  %12.6f  %10.6f\n", iteration, j, error_j);
    fprintf("         k  %12.6f  %10.6f\n", k, error_k);
    fprintf("         l  %12.6f  %10.6f\n",  l, error_l);
    fprintf("         m  %12.6f  %10.6f\n",  m, error_m);
    fprintf("         n  %12.6f  %10.6f  %10.6f\n", n, error_n, approx_error);
    fprintf("GAUSS-SEIDEL------------------------------------\n")
    
    iteration = iteration + 1;
end

% check
fprintf("\nCheck all equations equal zero with back substitution:\n");
fprintf("F1: %.6f\n", F1(j, k, l, m, n));
fprintf("F2: %.6f\n", F2(j, k, l, m, n));
fprintf("F3: %.6f\n", F3(j, k, l, m, n));
fprintf("F4: %.6f\n", F4(j, k, l, m, n));
fprintf("F5: %.6f\n", F5(j, k, l, m, n));






%%%%%% Jacobi
fprintf("\nJacobi Method\n")

% Initial guesses
j = 0;
k = 0;
l = 0;
m = 0;
n = 0;
error = 0.00001; % approximation error

% Initialize Variables
iteration = 1;
approx_error = Inf;

% Begin Iteration
fprintf("iter  unknown  value       error       maximum error \n");
while approx_error > error
    % Calculate new values using most recent value
    j_new = Fj(j, k, l, m, n);
    k_new = Fk(j, k, l, m, n);
    l_new = Fl(j, k, l, m, n);
    m_new = Fm(j, k, l, m, n);
    n_new = Fn(j, k, l, m, n);
    
    % Calculate error for each variable
    error_j = Ferror(j, j_new);
    error_k = Ferror(k, k_new);
    error_l = Ferror(l, l_new);
    error_m = Ferror(m, m_new);
    error_n = Ferror(n, n_new);

    % Continue iteration based on maximum error the set of variables
    approx_error = max([error_j, error_k, error_l, error_m, error_n]);

    j = j_new;    
    k = k_new;
    l = l_new;
    m = m_new;
    n = n_new;

    fprintf("%3i      j  %12.6f  %10.6f\n", iteration, j, error_j);
    fprintf("         k  %12.6f  %10.6f\n", k, error_k);
    fprintf("         l  %12.6f  %10.6f\n",  l, error_l);
    fprintf("         m  %12.6f  %10.6f\n",  m, error_m);
    fprintf("         n  %12.6f  %10.6f  %10.6f\n", n, error_n, approx_error);
    fprintf("JACOBI------------------------------------------\n")
    
    iteration = iteration + 1;
end

% check
fprintf("\nCheck all equations equal zero with back substitution:\n");
fprintf("F1: %.6f\n", F1(j, k, l, m, n));
fprintf("F2: %.6f\n", F2(j, k, l, m, n));
fprintf("F3: %.6f\n", F3(j, k, l, m, n));
fprintf("F4: %.6f\n", F4(j, k, l, m, n));
fprintf("F5: %.6f\n", F5(j, k, l, m, n));

