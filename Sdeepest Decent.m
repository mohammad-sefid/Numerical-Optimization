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

%% Steepest Descent Algorithm
g_k = gradient(f, X); % Gradient of f
func_eval = 0; % Number of function evaluations
grad_eval = 0; % Number of gradient evaluations
x(:, 1) = X0; % Start point
i = 1;

while true
    grad_f = subs(g_k, X, x(:, i)); % Gradient at x point
    grad_eval = grad_eval + 1;
    
    % Symbolic line search using Golden Section Search (GSS)
    alfa_opt(i) = GSS(f, X, x(:, i), grad_f); % Finding the optimum step size using GSS Method
    x(:, i + 1) = x(:, i) - alfa_opt(i) * grad_f; % Update the next point
    
    func_eval = func_eval + 1;
    
    epsilon = 10^(-3);
    if norm(x(:, i + 1) - x(:, i)) <= epsilon % Stop condition
        break
    end
    
    i = i + 1;
end

%% Display results
disp("Optimal solution: x* = ");
disp(double(x(:, end)));
disp("f(x*) = ");
disp(double(subs(f, X, x(:, end))));
disp("Number of function evaluations = " + num2str(func_eval));
disp("Number of gradient evaluations = " + num2str(grad_eval));
disp("Number of iterations = " + num2str(i - 1));

%% Golden Section Search (GSS) function
function gss_opt = GSS(f, X, x_i, grad_f)
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
        
        if subs(f, X, (x_i + a * grad_f)) < subs(f, X, (x_i + b * grad_f))
            XL = b;
        elseif subs(f, X, (x_i + a * grad_f)) > subs(f, X, (x_i + b * grad_f))
            XU = a;
        else
            XL = b;
            XU = a;
        end
    end
    
    gss_opt = (XL + XU) / 2;
end