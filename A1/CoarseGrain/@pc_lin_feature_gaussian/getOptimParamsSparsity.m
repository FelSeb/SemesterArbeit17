function [ this ] = getOptimParamsSparsity( this, suff_stats, data_provider, n_training_samples  )
% Computingoptimal theta_c and sigma^2 using a sparsity enforcing hyper
% prior. Optimization via EM. See Figueiredo.

% sigma2 is parametrized as follows:
% sigma2 = exp(l);

Max_iter = 5;
% sigma2_opt = this.Cov_matrix(1,1);
E_p_Q = zeros(1,Max_iter);

for iteration = 1:Max_iter
    sum_suff_stat2 = 0;
    sum_suff_stat1_feat = 0;
    sum_phi_phi = 0;
    for i = 1:n_training_samples
        [x_i, y_i] = data_provider.provideDataPoint(i);
        feat_fun_x_i = this.Feature_functions(x_i);
        
        sum_suff_stat2 = sum_suff_stat2 + suff_stats{2,i};
        sum_suff_stat1_feat = sum_suff_stat1_feat + suff_stats{1,i}' * feat_fun_x_i;
        sum_phi_phi = sum_phi_phi + feat_fun_x_i' * feat_fun_x_i;
        
    end
    sum_C = sum_suff_stat2...
        + 2 * sum_suff_stat1_feat * this.Coefficients...
        + this.Coefficients' * sum_phi_phi * this.Coefficients;
    
    if(iteration == 1)
        % number of features (= number of colums of feat_fun_x_i)
        [N_vars_coarse, M] = size(feat_fun_x_i);
    end
    
    % compute current optimal l
    l = log(sum_C/(M * n_training_samples));
    
    % compute current optimal sigma
    sigma2_opt = exp(l);
    % previousely:
    % sigma2_opt = -1/ (M * n_training_samples) * sum_C;

    % % exponential prior
    % gamma = 1;
    % V = diag(gamma * (1./this.Coefficients));
    
    % Jeffrey's prior (smarter computation according to Figueiredo 3.7)
    invcholV = diag(this.Coefficients);
    rhs = invcholV*(2/(n_training_samples) * sum_suff_stat1_feat');
    
    theta_opt = invcholV*(((invcholV*sum_phi_phi*invcholV/n_training_samples - eye(M)*sigma2_opt*2))\rhs);

%     % For testing purposes evaluate the lower bound: E_p(tau)[ Q(theata) ];
     V = diag((1./this.Coefficients).^2);
     E_p_Q(iteration) = -0.5*M*log(sigma2_opt) - 0.5/sigma2_opt*sum_C - 0.5*theta_opt'*V*theta_opt;
    
    this.Coefficients = theta_opt;
end

% For test purposes plot the lower bound over the iterations
plot(E_p_Q)

this.Cov_matrix = eye(N_vars_coarse) * sigma2_opt;
this.Cov_matrix_chol = eye(N_vars_coarse) * sqrt(sigma2_opt);
this.Det_cov_matrix = sigma2_opt^N_vars_coarse;
this.Log_det_cov_matrix = sum(2*log(diag(this.Cov_matrix_chol)));



end