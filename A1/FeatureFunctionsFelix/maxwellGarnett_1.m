function [feature] = maxwellGarnett_1( lambda, param_trafo , matrix_phase)
%Computes the effective conductivity using the Maxwell-Garnett Formula (see e.g. Torquato eq. 18.5)
%   lambda:         vector of conductivities
%   param_trafo:    transformation object

%ONLY FOR BINARY MATERIALS!

% Number of elements in the fine grid
n_lambda = length(lambda);

% low and high conductivity
lambda_hi = max(lambda);
lambda_lo = min(lambda);

% binary conductivity for the determination of volume fractions
lambda_binary = lambda > lambda_lo;

% Volume fraction of the high conducting phase
VolFrac_hi = sum(lambda_binary)/n_lambda;

if strcmp(matrix_phase, 'hi')
    %simplify notation
    lambdaMatrix = lambda_hi;
    lambdaInclusions = lambda_lo;
    inclusionVolFrac = 1 - VolFrac_hi;
    
    lambda_eff = lambdaMatrix*((lambdaMatrix + lambdaInclusions + inclusionVolFrac*(lambdaInclusions - lambdaMatrix))/...
        (lambdaMatrix + lambdaInclusions - inclusionVolFrac*(lambdaInclusions - lambdaMatrix)));
elseif strcmp(matrix_phase, 'lo')
    %simplify notation
    lambdaMatrix = lambda_lo;
    lambdaInclusions = lambda_hi;
    inclusionVolFrac = VolFrac_hi;
    
    lambda_eff = lambdaMatrix*((lambdaMatrix + lambdaInclusions + inclusionVolFrac*(lambdaInclusions - lambdaMatrix))/...
        (lambdaMatrix + lambdaInclusions - inclusionVolFrac*(lambdaInclusions - lambdaMatrix)));
else
    error('What is matrix phase for Maxwell-Garnett?')
end


% Transfromation to feature space
feature = param_trafo.getXfromElem2param(lambda_eff);

end

