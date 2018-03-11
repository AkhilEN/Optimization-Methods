function [ grad ] = grad_compute( fn , x )
% Computes the gradient of the function at x

% Step size for gradient calculation
step_size = 0.01;

% Length of variable vector
length_of_vector = length(x);

% Declaration of Zero Vector
e = zeros(length_of_vector,1);

% Pre - Initialization
grad = zeros(length_of_vector,1);

% Gradient calculation
for counter=1:length_of_vector
    e(counter)    = 1;
    grad(counter,1) = (fn(x + step_size*e) - fn(x - step_size*e))/(2*step_size);
    e(counter)    = 0;
end

end

