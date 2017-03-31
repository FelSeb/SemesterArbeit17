function [feature] = DEM_1( lambda, param_trafo , matrix_phase)
%Differential effective medium approximation, Torquato eq. 18.23


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

if strcmp(matrix_phase, 'lo')
    f = @(l) (lambda_hi - l)*sqrt(lambda_lo/l) - ...
        (1 - VolFrac_hi)*(lambda_hi - lambda_lo);
elseif strcmp(matrix_phase, 'hi')
    f = @(l) (lambda_lo - l)*sqrt(lambda_hi/l) - ...
        (1 - VolFrac_lo)*(lambda_lo - lambda_hi);
else
    error('DEM for high or low conducting phase as inclusion/matrix?')
end

lambda_eff = fzero(f, [lambda_lo lambda_hi]);

% Transfromation to feature space
feature = param_trafo.getXfromElem2param(lambda_eff);



end

