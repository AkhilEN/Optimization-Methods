function [ alpha ] = backtracking_line_search( fn, x, descent_direction )
% To calculate the step length in descent direction 

% Initial Step Length
alpha_init = 1;

% Parameter for Army Joe condition
c = 0.1;

% Correction of Step length
rho = 0.8;

% Descent Direction
descent_direction = descent_direction;

% Initialization
alpha = alpha_init;

% Gradient calculation
grad  = grad_compute(fn,x);

% Loop for calculating step length in desired descent direction
while (fn(x + alpha.*descent_direction) > fn(x) + c*alpha*(grad'*descent_direction))
    alpha = rho*alpha;
end


end

