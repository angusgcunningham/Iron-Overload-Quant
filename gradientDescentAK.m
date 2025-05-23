function [ropt, niter, gnorm, dr] = gradientDescentAK(r0, TE, a, y, mu, f_d2)

maxiter = 1000;         % Maximum number of allowed iterations
drmin = 1e-6;           % Minimum allowed perturbation
alpha = 1e-8;           % Step size (can be sensitive for different T2*)
tol = 1e-6;             % Termination tolerance

gnorm = inf;            % Onitialize gradient norm
r = r0;                 % Initilize result
niter = 0;              % Iteration counter
dr = inf;               % Perturbation

while and(gnorm >= tol, and(niter <= maxiter, dr >= drmin))
    gradientR = grad(TE, r, a, y, mu, f_d2);
    gnorm = norm(gradientR);
    rnew = r - alpha * gradientR;
    niter = niter + 1;
    dr = norm(rnew - r);
    r = rnew;
end
ropt = r;
niter = niter - 1;

end

function gradientR = grad(TE, r, a, y, mu, f_d2)
% %%%%%% INSERT YOUR CODE HERE % %%%%%%
    % 1) Residual: (s - a exp(-r TE))
    TE = TE(:);
    residual = y - a * exp(-r * TE);     
    
    % 2) dataTerm: a [TE .* exp(-r TE)]^T residual
    dataTerm = ( a * (TE .* exp(-r * TE)) )' * residual; 
    regTerm = mu * (r - f_d2);
    
    % 4) Final gradient
    gradientR = dataTerm + regTerm;
end
