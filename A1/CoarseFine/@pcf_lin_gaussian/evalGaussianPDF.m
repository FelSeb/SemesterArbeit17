function [val, Y, dval] = evalGaussianPDF(this, y, X, eval_log, varargin)
% PCF
% X must be in vector format

if(numel(varargin)>0)
    compute_log_grad = varargin{1};
else
    compute_log_grad = 0;
end

k = length(y);

% Evaluate coarse model
[Y, this.Solver] = this.Solver.evaluateModel(X,1);

% % Evaluate coarse model and return gradient of the solution vector
% [Y2,~,dYdX] = this.Solver.evaluateModel(X,1,1);

% This guarantees postivity of the exponent
% Note: a'inv(S)a = a'inv(R'R)a = a'inv(R)inv(R')a = Rit_a' * Rit_a
vec = (y - this.Offset - this.Interpolation_matrix * Y);
Rit_a = this.S_chol'\vec;

if(~eval_log)
    % if eval_log = False evalueate the PDF
    %     mean = Offset + Interpolation_matrix * Y;
    %     val = (2*pi)^(-k/2) * det(Cov_matrix)^-0.5 * exp(-0.5* (y - mean)' * inv(Cov_matrix) * (y - mean));
    val = (2*pi)^(-k/2) * this.Det_S^-0.5 * exp(-0.5 * (Rit_a' * Rit_a) );
else
    % if eval_log = True, evalueate the log PDF
    val = (-k/2)*log(2*pi) -0.5* this.Log_det_S -0.5 * (Rit_a' * Rit_a);
    
%     % Computation of the gradient of the log PDF w.r.t. material parameters
%     % Bad idea (may be used for testing ...compare result to method below)
%     if(compute_log_grad)
%         dval2 = vec' * (this.S_chol\ (this.S_chol'\ (this.Interpolation_matrix * dYdX)));
%     end
    
    % Computation gradient of the log PDF w.r.t. material parameters
    if(compute_log_grad)
        dval = getPCFGrad(this, Y, y);
    end
    
end

end


function dval = getPCFGrad(pcf, u, y)
% All equations corresponding to nondir nodes
nondir = pcf.Coarse_grid.Node2eq(pcf.Coarse_grid.All_nondir_nodes);

% Computation of lambda_prefac from the adjoint method (efficient method
% for the computation of the gradient of pcf).
WtSinv = (pcf.S_chol\(pcf.S_chol'\pcf.Interpolation_matrix(:,nondir)))';
%WtSinv = (S\W)';
lambda_prefac = pcf.Solver.C'\ (WtSinv * (y - pcf.Offset - pcf.Interpolation_matrix(:,nondir)*u(nondir)) );

Nx = numel(pcf.Solver.C_grad);
syssize = size(pcf.Solver.C_grad{1},1);

C_gradu = zeros(syssize,Nx);

% Should be fast since each C_grad{i} has only 8 non-zero entries
for ix = 1:Nx
    C_gradu(:,ix) = pcf.Solver.C_grad{ix} * u(nondir);
end

dval = - lambda_prefac' * (C_gradu - pcf.Solver.Rhs_grad);

end







