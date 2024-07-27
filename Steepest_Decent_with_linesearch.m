clc; clear all; close all;

% Define symbolic variables
syms x1 x2 x3 x4 alpha

% Define function 1
f = 100 * (x2 - (x1)^2)^2 + (1 - x1)^2;
x = [1 2];
X = [x1, x2];

% Define function 2 
% f = (x1 + 10*x2)^2 + 5*(x3 - x4)^2 + (x2 - 2*x3)^4 + 10*(x1 - x4)^4;
% x = [1 2 2 2];
% X = [x1, x2, x3, x4];

% Set tolerance and initialize counters
epsilon = 10^-3;
g = gradient(f);  % Gradient vector
fun_eval = 0;  % Number of function evaluations
grad_eval = 0;  % Number of gradient evaluations
flag = 1;
iterations = 0;  % Number of iterations

% Optimization loop
while flag == 1
    iterations = iterations + 1;
    grad_i = double(subs(g, X, x)).';  % Gradient vector at current point
    grad_eval = grad_eval + 1;
    p_i = -grad_i;  % Search direction
    f_i = subs(f, X, x + alpha * p_i);  % f(x + alpha*p)
    [alpha_i, f_e, g_e] = linesearch(f_i, 2);  % Line search
    fun_eval = fun_eval + f_e;  % Update function evaluation count
    grad_eval = grad_eval + g_e;  % Update gradient evaluation count
    x_prev = x;  % Previous x
    x = x + alpha_i * p_i;  % Update x
    if double(norm(x - x_prev)) <= epsilon  % Stop condition
        flag = 0;
        break;
    end
end

% Display results
disp("x* = ");
disp(double(x));
disp("f(x*) = ");
disp(double(subs(f, X, x)));
disp("Number of function evaluations = " + num2str(fun_eval));
disp("Number of gradient evaluations = " + num2str(grad_eval));
disp("Number of iterations = " + num2str(iterations));