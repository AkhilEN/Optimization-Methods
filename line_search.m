function [ min ] = line_search( fn , xguess )
% Calculates Local minimum using line search methods

% Iteration number
k = 0;

% Maximum number of iterations
N = 1000;

% Tolerence on square of gradient norm
epsilon = 10^-4;

% Gradient calculation
grad = grad_compute(fn , xguess);

% Assigning initial guess
x = xguess;

% Loop to calculate the local minima
while ((norm(grad)>epsilon) && (k<=N))
    p = -grad_compute(fn , x);
    alpha = backtracking_line_search(fn , x , p);
    x = x + alpha.*p;
    k = k +1;
    grad = grad_compute(fn , x);
end

min = x;


end

