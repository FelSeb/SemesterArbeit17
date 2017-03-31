function [ this ] = getOptimSDiag( this, suff_stats, data_provider, n_training_samples  )

optim_S_diag = 0;

for i = 1:n_training_samples
    [x_i, y_i] = data_provider.provideDataPoint(i);
    % Compute W * <Y(X_i)>
    WY = this.Interpolation_matrix * suff_stats{3,i};
    
    % Compute (y_i' - mu')
    y_i_min_mu = y_i - this.Offset;
    
    % Compute (y_i' - mu' - <Y(X_i)>'*W')
    term1 = y_i_min_mu' - WY';

    
    WYYWt = zeros(length(y_i),1);
    for ii = 1:length(y_i)
        WYYWt(ii) = this.Interpolation_matrix(ii,:) * suff_stats{4,i} * this.Interpolation_matrix(ii,:)';
    end
    
    % compose everything
    optim_S_diag = optim_S_diag ...
        + y_i_min_mu .* term1'...
        - WY .* y_i_min_mu...
        + WYYWt;
    
end

optim_S_diag = optim_S_diag / n_training_samples;

this.S = diag(optim_S_diag);
this.S_chol = diag(sqrt(optim_S_diag));
this.Det_S = prod(optim_S_diag);
this.Log_det_S = sum(log(optim_S_diag));

end

