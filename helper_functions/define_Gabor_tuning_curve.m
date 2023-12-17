% run this script to put anonymous functions for a Gabor tuning curve and
% it's first derivative into the workspace

% Define function for Gabor tuning curve
gabor_tuning_curve = @(x,p) p(1) + p(2)*exp( -(x-p(3)).^2 / (2*p(4)^2) ) .* cos(2*pi*p(5)*(x-p(3))+p(6));

% Define function for derivative of Gabor tuning curve
first_deriv_of_GTC = @(x,p) -p(2)*exp(-0.5*((x - p(3))/p(4)).^2 ).* ...
    ( (x - p(3))/(p(4)^2).*cos(2*pi*p(5)*(x - p(3)) + p(6))...
    + (2*pi*p(5))*sin(2*pi*p(5)*(x - p(3)) + p(6)) );