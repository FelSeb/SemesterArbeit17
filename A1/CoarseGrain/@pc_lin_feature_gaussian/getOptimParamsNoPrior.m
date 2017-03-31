function [ this ] = getOptimParamsNoPrior( this, suff_stats, data_provider, n_training_samples )
% PC

warning('sigma2 is not expressed as exp(log(sigma2)) here. this may lead to an negative optimal sigma2. This should be changed')

%N_elem_coarse = this.Coarse_grid.N_elem;
N_vars_coarse = numel(suff_stats{1,1});

% Precomputation of certain values. In a later version of the code the
% training data should be loaded in batches and not all at once.
sum_Phi_Phi = 0;
sum_suff_stat2 = 0;
sum_feat_suff_stat1 = 0;
for i = 1:n_training_samples
    [x_i, y_i] = data_provider.provideDataPoint(i);
    feat_fun_x_i = this.Feature_functions(x_i);
    sum_Phi_Phi = sum_Phi_Phi +  feat_fun_x_i' * feat_fun_x_i;
    
    % compute sums over sufficient statistics
    sum_suff_stat2 = sum_suff_stat2 + suff_stats{2,i};
    
    sum_feat_suff_stat1 = sum_feat_suff_stat1 + feat_fun_x_i' * suff_stats{1,i};
    
end

% Using no piror we obtain the remaining parameters theta_c via the
% solution of a linear system of equations:
R = chol(sum_Phi_Phi);
this.Coefficients = R\(R'\sum_feat_suff_stat1);


sum_coeff_feat_suff_stat1 = this.Coefficients' * sum_feat_suffstat1;

% Computation of variance parameter
sigma2 = 1/(n_training_samples*N_vars_coarse) * ...
    (sum_suff_stat2  - 2 * sum_coeff_feat_suff_stat1 + this.Coefficients' * sum_Phi_Phi * this.Coefficients );


this.Cov_matrix = eye(N_vars_coarse) * sigma2;
this.Cov_matrix_chol = eye(N_vars_coarse) * sqrt(sigma2);
this.Det_cov_matrix = sigma2^N_vars_coarse;
this.Log_det_cov_matrix = sum(2*log(diag(this.Cov_matrix_chol)));



% Employing a prior depending on theta_c we must compose a system of linear
% equations:
% TO DO...

end






