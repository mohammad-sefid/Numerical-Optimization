clc
clear
close all

syms x1 x2 x3 x4 alfa;

%% Asking the user to choose a function
fprintf('Choose a function:\n');
fprintf('1. f = 100*((x2-x1^2)^2)+(1-x1)^2 X0=[1;2]; chong 5.1\n');
fprintf('2. f = (x1+10*x2)^2+5*(x3-x4)^2+(x2-2*x3)^4+10*(x1-x4)^4 X0=[1;2;2;2]; chong 9.1\n');

choice = input('Enter a function (1 or 2): ');

%% Choose function
switch choice
    case 1 
        f = 100*((x2-x1^2)^2)+(1-x1)^2;
        X0 = [1; 2];
        X = [x1; x2];
    case 2
        f = (x1+10*x2)^2+5*(x3-x4)^2+(x2-2*x3)^4+10*(x1-x4)^4;
        X0 = [1; 2; 2; 2];
        X = [x1; x2; x3; x4];
end

%% Newton Algorithm
grad_func = gradient(f, X); % Gradient of function
hessian_func = hessian(f, X); % Hessian of function
func_eval = 0; % Number of function evaluations
grad_eval = 0; % Number of gradient evaluations
hessian_eval = 0; % Number of hessian evaluations
x(:, 1) = X0; % Start point
i = 1;

while true
    gradient = subs(grad_func, X, x(:, i)); % Gradient at x point
    hessian = subs(hessian_func, X, x(:, i)); % Hessian at x point
    gradient_evaluation = grad_eval + 1;
    hessian_evaluation = hessian_eval + 1;
    
    % Check if Hessian is not Positive Definite (P.D)
    if min(double(eig(hessian))) <= 0
        lambda = abs(min(double(eig(hessian))));
        hessian = hessian + (lambda + 10^(-6)) * eye(size(hessian)); % Adjust hessian to be P.D.
    end
    
    hessian_inv = inv(hessian); % Inverse of Hessian matrix
    p_k = -hessian_inv * gradient; % Calculate the search direction vector p
    
    % Symbolic line search using Golden Section Search (GSS)
    alfa_opt = GSS(f, X, x(:, i), p_k);
    
    % Update x
    x(:, i + 1) = x(:, i) + alfa_opt * p_k;
    
    % Check convergence
    epsilon = 10^(-3);
    if norm(x(:, i + 1) - x(:, i)) <= epsilon
        break
    end
    
    i = i + 1;
end

%% Display results
disp("Optimal solution: x* = ");
disp(double(x(:, end)));
disp("Objective function value: f(x*) = ");
disp(double(subs(f, X, x(:, end))));
disp("Number of function evaluations: " + num2str(func_eval));
disp("Number of gradient evaluations: " + num2str(gradient_evaluation));
disp("Number of hessian evaluations: " + num2str(hessian_evaluation));
disp("Number of iterations: " + num2str(i - 1));

%% Golden Section Search (GSS) function
function gss_opt = GSS(f, X, x_i, p_k)
    syms alfa;
    XL = 0;
    XU = 2;
    ro = 0.382;
    acc = 10^(-5);
    N = ceil(log(acc / (XU - XL)) / log(1 - ro)); % Number of iterations
    
    % Golden Section Search (GSS) algorithm
    for j = 1:N
        d = (1 - ro) * (XU - XL);
        a = XL + d;
        b = XU - d;
        
        if subs(f, X, (x_i + a * p_k)) < subs(f, X, (x_i + b * p_k))
            XL = b;
        elseif subs(f, X, (x_i + a * p_k)) > subs(f, X, (x_i + b * p_k))
            XU = a;
        else
            XL = b;
            XU = a;
        end
    end
    
    gss_opt = (XL + XU) / 2;
end
