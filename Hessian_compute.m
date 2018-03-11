function [ Hessian ] = Hessian_compute( fn , x )
% Computes the Hessian of the function fn at x

% Step size for gradient calculation
step_size = 0.01;

% Length of variable vector
length_of_vector = length(x);

% Declaration of Zero Vectors
ei = zeros(length_of_vector,1);
ej = zeros(length_of_vector,1);

% Hessian calculation
for counter1=1:length_of_vector
    ei(counter1) = 1;
    for counter2=1:length_of_vector
        ei(counter1) = 1;
        ej(counter2) = 1;
        Hessian(counter1,counter2) = (fn( x + step_size*ei + step_size*ej ) - fn ( x + step_size*ei ) - fn( x + step_size*ej) + fn( x ))/(step_size^2);
        ej(counter2) = 0;
    end
    ei(counter1) = 0;
end

end

