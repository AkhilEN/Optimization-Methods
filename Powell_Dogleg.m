function [ p ] = Powell_Dogleg( fn,x,Trustradius )
% Calculates the Step Length required in the Trust Region method as per the
% Powell Dogleg method

% Calculation of Newton Step
    Gradient = grad_compute(fn,x);
    Hessian  = Hessian_compute(fn,x);
    pNewton  = - Hessian\Gradient;
    % Newton Step Check (1st Priority)
    if norm(pNewton)<= Trustradius
        % Newton Step below Trustradius
        p = pNewton;
        return;
    end
% Calculation of Cauchy Step
    % Test for Determining positive Definiteness
    Test = Gradient'*Hessian*Gradient; 
    % Assignment of Cauchy Step
    if Test>0
        pCauchy = - ((Gradient'*Gradient)/Test)*Gradient;
    else
        % Steepest Descent
        pCauchy = - Trustradius/norm(Gradient) * Gradient;
    end
    % Cauchy Step Check (2nd Priority)
    if norm(pCauchy)>=Trustradius
        p = Trustradius/norm(pCauchy) * pCauchy ;
        return;
    end
% Integrated Step Calculation (Under ||p|| = Trustradius)
    eta = (-(pNewton-pCauchy)'*pCauchy...
    + ((((pNewton - pCauchy)'*pCauchy)^2)...
    + (Trustradius^2 - norm(pCauchy)^2)*(norm(pNewton-pCauchy))^2)^0.5)...
    /norm(pNewton-pCauchy)^2 ;
    % Integrated Step Check
    if eta>0
        p = eta*pNewton + (1-eta)*pCauchy;
        return;
    end
    
    

end

