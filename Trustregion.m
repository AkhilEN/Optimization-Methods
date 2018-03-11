function [ fnmin,xmin ] = Trustregion( fn , x0 )
% Finds the minima of the function fn using trust region methods 
% Input Arguments are the function under consideration(fn) and the starting
% point(x0)

% Maximum Step Length
Delta0 = 10;
% All Step Lengths lie between 0 and Delta0

% Initial Trustradius
Trustradius = Delta0/2;

% Defining the minimum reduction ratio for us to consider the decrease
% significant
% If less than this we do not move from current point and decrease the
% Trust Radius
% Generally from  0 to 0.25
eta = 0.25;

% Tolerence Limit on Convergence
epsilon = 10^-05;

% Initial Tolerence Limit (For starting the while loop)
tol = 1;

% Defining Zero Vector
% Number of Variables
    len  = length(x0);
    null = zeros(len,1);

% Initial Assignment
x = x0;

while abs(tol)>epsilon
    % Construction of the Quadratic Approximation of the function around
    % the current point x
    m = @(p) fn(x)+((grad_compute(fn,x))'*p)...
    +(0.5*(p'*Hessian_compute(fn,x)*p));
    % Calculation of Step Length using Powell-Dogleg Method
    p = Powell_Dogleg(fn,x,Trustradius);
    % Calculation of Reduction Ratio    
    rho = (fn(x) - fn(x+p))/(m(null) - m(p));
    % Checking Based on Reduction Ratio
    if rho < 0.25
        Trustradius = 0.25 * norm(p);
    else
        if rho > 0.75 && norm(p) == Trustradius
            % Function Approximation is good and we are achieving good
            % reduction. Hence increasing Trust Region
            Trustradius = min(2*Trustradius,Delta0);
        else
            % The Trustradius is just fine
            Trustradius = Trustradius;
        end
    end
    % Checking for sufficient Decrease of function
    if rho > eta
        % Decrease achieved greater than set minimum
        x = x + p;
        % Calculation on tolerence for convergence
        tol = fn(x) - fn(x-p);
    else
        % Sufficient Decrease not achieved
        x = x;
        % No Tolerence Calculation as no step has been achieved
    end    
end
% Return the minima point
xmin = x;
% Return the minimum function value
fnmin = fn(x); 

end

