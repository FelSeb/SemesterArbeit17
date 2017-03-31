function [ this ] = getOptimSFull( this, suff_stats, data_provider, n_training_samples  )
warning('this functions needs to be checked extensively')





optim_S = 0;

for i = 1:n_training_samples
    [x_i, y_i] = data_provider.provideDataPoint(i);
    % Compute W * <Y(X_i)>
    WY = this.Interp_matrix * suff_stats{3,i};
    
    % Compute (y_i' - mu')
    y_i_min_mu = yi - this.Offset;
    
    % Compute (y_i' - mu' - <Y(X_i)>'*W')
    term1 = y_i_min_mu' - WY';

    % compose everything
    optim_S = optim_S ...
        + y_i_min_mu * term1...
        - WY * y_i_min_mu'...
        + this.this.Interp_matrix * suff_stats{4,i} * this.this.Interp_matrix';
    
end

optim_S = optim_S / n_training_samples;

this.S = optim_S;
this.S_chol = chol(optim_S);
this.Det_S = prod(diag(this.S_chol))^2;
this.Log_det_S = 2*sum(log(diag(this.S_chol)));


end

