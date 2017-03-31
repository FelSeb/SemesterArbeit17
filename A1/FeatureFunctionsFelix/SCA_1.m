function [ feature ] = SCA_1( lambda, param_trafo )
% It is assumed that the paramter field takes one of two values.

% Number of elements in the fine grid
n_lambda = length(lambda);

% low and high conductivity
lambda_hi = max(lambda);
lambda_lo = min(lambda);

% binary conductivity for the determination of volume fractions
lambda_binary = lambda > lambda_lo;

% Volume fractions
VolFrac_hi = sum(lambda_binary)/n_lambda;
VolFrac_lo = 1-VolFrac_hi;


alpha = lambda_lo*(2*VolFrac_lo - 1) + lambda_hi*(2*VolFrac_hi - 1);

lambda_eff = 0.5*(alpha + sqrt(alpha^2 + 4*lambda_lo*lambda_hi));


% Transfromation to feature space
feature = param_trafo.getXfromElem2param(lambda_eff);

end

