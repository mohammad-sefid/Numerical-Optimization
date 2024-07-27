function CA5_BFGS_40206864()
    clc; clear all; close all;

    % Define the objective function
    syms x1 x2
    f_sym = 100*(x2 - (x1)^2)^2 + (1 - x1)^2; % Rosenbrock function

    % Convert symbolic function to MATLAB function handle
    f = matlabFunction(f_sym, 'Vars', [x1, x2]);
    grad_f_sym = gradient(f_sym, [x1, x2]);
    grad_f = matlabFunction(grad_f_sym, 'Vars', [x1, x2]);

    % Set the initial point
    x0 = [-1.2; 1];  % Common starting point for Rosenbrock function

    % Define the termination criteria
    tol = 1e-6;
    max_iter = 1000;

    % Initialize the BFGS algorithm
    x = x0;
    B = eye(length(x0)); % Initial Hessian approximation
    k = 0;
    f_val = f(x(1), x(2));
    grad_val = grad_f(x(1), x(2));
    func_evals = 1;
    grad_evals = 1;

    % BFGS algorithm
    while norm(grad_val) > tol && k < max_iter
        % Compute the search direction
        p = -B \ grad_val;

        % Perform line search
        [alpha, f_new, grad_new, fe, ge] = linesearch(f, grad_f, x, p);
        func_evals = func_evals + fe;
        grad_evals = grad_evals + ge;

        % Update the solution
        x_new = x + alpha * p;

        % Update the Hessian approximation
        s = x_new - x;
        y = grad_new - grad_val;
        rho = 1 / (y' * s);
        B = (eye(length(x0)) - rho * (s * y')) * B * (eye(length(x0)) - rho * (y * s')) + rho * (s * s');

        % Update the iteration
        x = x_new;
        f_val = f_new;
        grad_val = grad_new;
        k = k + 1;
    end

    % Display the results
    fprintf('Optimal solution: x* = [%.4f, %.4f]\n', x(1), x(2));
    fprintf('Objective function value: f(x*) = %.4f\n', f_val);
    fprintf('Number of function evaluations: %d\n', func_evals);
    fprintf('Number of gradient evaluations: %d\n', grad_evals);
    fprintf('Number of iterations: %d\n', k);
end

function [alpha_star, f_new, grad_new, f_eval, g_eval] = linesearch(f, grad_f, x, p)
    % Parameters for line search
    c1 = 0.1;
    c2 = 0.9;
    alpha = 1.0;  % Initial step length
    alpha_low = 0;
    alpha_high = inf;
    max_iter = 100;
    iter = 0;

    % Evaluate the objective function and its gradient at the initial point
    f_x = f(x(1), x(2));
    grad_x = grad_f(x(1), x(2));

    f_eval = 1;
    g_eval = 1;

    % Perform Wolfe condition line search
    while iter < max_iter
        % Compute the new candidate point and function value
        x_new = x + alpha * p;
        f_new = f(x_new(1), x_new(2));
        grad_new = grad_f(x_new(1), x_new(2));

        % Check the Wolfe conditions
        if f_new > f_x + c1 * alpha * (grad_x' * p)
            alpha_high = alpha;
        elseif grad_new' * p < c2 * grad_x' * p
            alpha_low = alpha;
        else
            alpha_star = alpha;
            return;
        end

        if alpha_high < inf
            alpha = (alpha_low + alpha_high) / 2;
        else
            alpha = 2 * alpha;
        end

        iter = iter + 1;
        f_eval = f_eval + 1;
        g_eval = g_eval + 1;
    end

    alpha_star = alpha;
end