clc
clear
close all

%% Asking the user to choose a function
fprintf('Choose a function:\n');
fprintf('1. f = 100*((x2-x1^2)^2)+(1-x1)^2 X0=[1;2]; chong 5.1\n');
fprintf('2. f = (x1+10*x2)^2+5*(x3-x4)^2+(x2-2*x3)^4+10*(x1-x4)^4 X0=[1;2;2;2]; chong 9.1\n');

choice = input('Enter a function (1 or 2): ');

%% Choose function
switch choice
    case 1 
        f = @(X) 100*((X(2)-X(1)^2)^2)+(1-X(1))^2;
        X0 = [1; 2];
    case 2
        f = @(X) (X(1)+10*X(2))^2 + 5*(X(3)-X(4))^2 + (X(2)-2*X(3))^4 + 10*(X(1)-X(4))^4;
        X0 = [1; 2; 2; 2];
end

%% Steepest Descent Algorithm
func_eval = 0; % Number of function evaluations
grad_eval = 0; % Number of gradient evaluations
x(:, 1) = X0; % Start point
i = 1;

while true
    % Compute gradient numerically using finite difference
    grad_f = numerical_gradient(f, x(:, i)); % Gradient at x point
    grad_eval = grad_eval + 1;
    
    % Symbolic line search using Golden Section Search (GSS)
    alfa_opt(i) = GSS(f, x(:, i), grad_f); % Finding the optimum step size using GSS Method
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
disp(x(:, end));
disp("f(x*) = ");
disp(f(x(:, end))); % Evaluating function at optimal point
disp("Number of function evaluations = " + num2str(func_eval));
disp("Number of gradient evaluations = " + num2str(grad_eval));
disp("Number of iterations = " + num2str(i - 1));

%% Golden Section Search (GSS) function
function gss_opt = GSS(f, x_i, grad_f)
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
        
        % Evaluate function at a and b
        fa = f(x_i + a * grad_f);
        fb = f(x_i + b * grad_f);
        
        if fa < fb
            XL = b;
        elseif fa > fb
            XU = a;
        else
            XL = b;
            XU = a;
        end
    end
    
    gss_opt = (XL + XU) / 2;
end

%% Numerical gradient calculation using finite difference
function grad_f = numerical_gradient(f, x)
    h = 1e-5; % Small perturbation
    n = length(x);
    grad_f = zeros(n, 1);
    
    for i = 1:n
        x1 = x;
        x1(i) = x1(i) + h;
        grad_f(i) = (f(x1) - f(x)) / h;
    end
end
