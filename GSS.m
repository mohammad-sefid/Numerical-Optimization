clc 
clear 
close all

% input parameters
disp("Function is defined as: A*x^4 + B*x^3 + C*x^2 + D*x + E"); 
A = input("enter A \n"); 
B = input("enter B \n"); 
C = input("enter C \n"); 
D = input("enter D \n"); 
E = input("enter E \n"); 
disp("range is: (L0,U0)"); 
L0 = input("enter L0 \n");
U0 = input("enter U0 \n");
accuracy = input("enter accuracy \n");

% define function
syms x
f = A*x^4+B*x^3+C*x^2+D*x+E; 

% finding N: number of iterations needed
p = 0.382; % ro in the GSS is constant
flag = 1;
N=1;

% computing number of iterations needed to get desired accuracy
while(flag == 1)
    n_accuracy = (U0-L0)*((1-p)^N); % n_accuracy is accuracy with current number of iterations
    if n_accuracy <= accuracy
       break % stop when the accuracy is enough
    end
    N = N+1; % increase N when the accuracy isn't enough
end

% Algorithm iterations
for i =1:N
    delta = (U0-L0)*p; % distance between ends of range (L0 , U0) and selected points (L1 , U1) 
    L1 = L0 + delta;% choosing  L1 point
    U1 = U0 - delta;% choosing  U1 point
    f_L1 = double(subs(f,L1)); %  caluculate f(L1)
    f_U1 = double(subs(f,U1));%  caluculate f(U1)
    if f_L1 > f_U1
       L0 = L1; % new range is (L1,U0), as in the Algorithm process
    else
       U0 = U1; % new range is (L0,U1), as in the Algorithm process
    end
end
x_star = (L0+U0)/2; % estimated answer

% displaying results
disp("answer is x* and a0 <= x* <= b0 where:");
disp("L0 = "+num2str(L0));
disp("U0 = "+num2str(U0));
disp("estimated answer is : "+num2str(x_star));
disp("f(L0) = "+num2str(double(subs(f,L0))));
disp("f(U0) = "+num2str(double(subs(f,U0))));
disp("f(x*) = "+num2str(double(subs(f,x_star))));