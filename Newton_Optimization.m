function [ x_optim ] = Newton_Optimization( fn , x0)
% Calculates the minimum value of the function fn using Newton approach

% Initial Guess - x0
xold = x0;

% Tolerence Limits
epsilon1 = 10^-3;
epsilon2 = 10^-3;

% Evaluating gradient and Hessian at initial guess
grad    = grad_compute(fn , x0);
Hessian = Hessian_compute(fn , x0);

% Iteration Number
k = 1;

p = grad;

while norm(p) > epsilon1 && norm(grad) > epsilon2
    
    grad    = grad_compute(fn , xold);
    Hessian = Hessian_compute(fn , xold);
    if (rcond(Hessian) > 1.5 || rcond(Hessian) < 0.5)
        break;
    else
        p    = Hessian\(-1*grad);
        xnew = xold + p;
    end
    xold = xnew;
    k    = k + 1 ;
end

x_optim = xold;

disp(k);
end

