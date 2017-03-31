% This is a file contining information on the different data types occuring
% in the code:

% Elem2param can take one of the three forms (lambdaELEMENT_PARAMNUMBER):
% 1.)   [lambda1_1, lambda2_2, lambda3_3, lambda4_4,... ]
%
% 2.)   [lambda1_1, lambda1_2; lambda2_1, lambda2_2; lambda3_1,... ]'
%
% 3.)   [[lambda1_1, lambda1_2; lambda1_2, lambda1_3],
%        [lambda2_1, lambda2_2; lambda2_2, lambda2_3], ...]


% The corresponding parametrization of lambda through X is done as follows:
% 1.)   lambda_1 = exp(2 * X_1)
% 2.)   lambda_1 = exp(2 * X_1)
%       lambda_2 = exp(2 * X_2)
% 3.)   lambda_1 = exp(2 * X_1)
%       lambda_2 = exp(2 * X_2) + X_3^2
%       lambda_3 = exp(2 * X_1) * X_3 
        

