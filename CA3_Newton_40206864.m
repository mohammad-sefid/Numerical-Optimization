clc
clear
close all

syms x1 x2 x3 x4 alfa;
%% asking the user to choose a function
fprintf('Choose a function:\n');
fprintf('1. f = 100*((x2-x1^2)^2)+(1-x1)^2 X0=[1;2]; chong 5.1\n');
fprintf('2. f=(x1+10*x2)^2+5*(x3-x4)^2+(x2-2*x3)^4+10*(x1-x4)^4 X0=[1;2;2;2]; chong 9.1\n');

choice = input('Enter a function (1 or 2):');

%% choose function
switch choice
    case 1 
        f = 100*((x2-x1^2)^2)+(1-x1)^2;
        X0=[1;2];
        X = [x1;x2];
    case 2
        f=(x1+10*x2)^2+5*(x3-x4)^2+(x2-2*x3)^4+10*(x1-x4)^4;
        X0=[1;2;2;2];
        X=[x1;x2;x3;x4];
end
%% Newton Algorithm
grad_func=gradient(f,X); %gradient of function
hessian_func=hessian(f,X); %hessian of function
func_eval=0;  % number of function evaluations
grad_eval=0;  % number of gradient evaluations
hessian_eval=0; % number of hessian evaluations
x(:,1)= X0 ;  % start point
i = 1;
while 1
    gradient=subs(grad_func,X,x(:,i)); % gradient at x point
    hessian=subs(hessian_func,X, x(:,i)); % hessian at x point
    gradient_evaluation=grad_eval+1;
    hessian_evaluation=hessian_eval+1;
    if hessian<=0 % checking hessian matrix is not Positive Definite(P.D)
        lambda=abs(eig(hessian));
        hessian=hessian+(lambda+10^(-6)); % adjust the hessian matrix to be Positive Definite(P.D)
    end
    hessian_inv=inv(hessian); % Inverse of hessian matrix
    p_k=-hessian_inv*gradient; % calculate the vector p
    alpha_k=subs(f,X,(x(:,i)+alfa*p_k));
    func_eval=func_eval+1;
    alfa_opt=GSSMethod(alpha_k,0,2); % finding the optimum step size using GSS Method
    x(:,i+1)=x(:,i)+alfa_opt*p_k; % the next point
    epsilon=10^(-3);
    if norm(x(:,i+1)-x(:,i))<=epsilon % stop condition
        break
    end
    i=i+1;
end
%% Display results
disp("x* = ");
disp(double(x(:,end)));
disp("f(x*) = ");
disp(double(subs(f,X,x(:,end))));
disp("number of function evaluations = "+num2str(func_eval));
disp("number of gradient evaluations = "+num2str(gradient_evaluation));
disp("number of hessian evaluations = "+num2str(hessian_evaluation));
disp("number of iterations = "+num2str(i-1));
