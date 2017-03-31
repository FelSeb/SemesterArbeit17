function [val, dlogval] = evalGaussianPDF(this, X, lambda, eval_log, varargin)
% PC
% X and x_tensor are the coarse and fine inputs.
%
% X_tensor must be in vector format
% x_tensor must have one of the Elem2param formats. ???

if(numel(varargin)>0)
    compute_log_grad = varargin{1};
else
    compute_log_grad = 0;
end

k = length(X);

% convert X to column vector
if(length(X) == numel(X))
    X = reshape(X,[numel(X),1]);
else
    error('X has wrong format');
end
% Note: a'inv(S)a = a'inv(R'R)a = a'inv(R)inv(R')a = Rit_a' * Rit_a
Rit_a = this.Cov_matrix_chol'\(X - this.Feature_functions(lambda) * this.Coefficients);


if(~eval_log)
    % if eval_log = False evalueate the PDF
    %     mean = weights' * Feature_functions(x);
    %     Cov_matrix = sigma^2 * I;
    %     val = 1/sqrt(2*pi * det(Cov_matrix)) * exp(-0.5* (x - mean)' * inv(Cov_matrix) * (x - mean));
    val = (2*pi)^(-k/2) * this.Det_cov_matrix^-0.5 * exp(-0.5 * (Rit_a' * Rit_a) );
    dlogval = [];
else
    % if eval_log = True evalueate the log PDF
    val = (-k/2)*log(2*pi) -0.5* this.Log_det_cov_matrix - 0.5 * (Rit_a' * Rit_a);
    
    % Computation of the gradient of the log PDF w.r.t. material parameters
    if(compute_log_grad)
        dlogval = -(this.Cov_matrix_chol\(Rit_a))';
    else
        dlogval = [];
    end
end

end