function [ SuffStats ] = getSuffStats( this )
% Computing a MCMC Method to compute the sufficient statistics.

n_coarse_elem = this.PCF.Coarse_grid.N_elem;

E2Pdim = size(this.PCF.Coarse_grid.Elem2param);

% This for loop may be parallelized!!!
for i = 1:this.N_training_samples
    % set training data x and y of current batch
    [x_i, y_i] = this.Data_provider.provideDataPoint(i);
    
    % Target distribuiton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    warning('shouldnt this be: this.logQopt_i(X,x_i,y_i, 1)')
    LOGPDF = @(X) logQopt_i(this, X, x_i, y_i, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set initial guess. Using the mean of pc seems to be a convenient
    % choice.
    init = (this.PC.Feature_functions(x_i) * this.PC.Coefficients)';
    
    % Determine appropriate stepwidth
    this.Sampler = this.Sampler.optimizeStepwidth(init, LOGPDF);

    % Generate all X-samples
    tic;
    [smpl,accept] = this.Sampler.generateSamples(LOGPDF, init);
    toc;
    
    % Derive <X>, <|X|^2>, <Y(X)>, <Y(X)Y(X)'>
    %all_X = smpl;
    all_X = smpl.samples;
    all_Y = smpl.info;
    
    %%%%%%%%%%%%%%%%%%% for test purposes &%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dvals = smpl.info2;
%     auto_cov = this.Sampler.estimateAutoCovariance(all_X,1/5);
%     accept
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SuffStats{1,i} = mean(all_X',2);
    SuffStats{2,i} = 0;
    for j = 1: length(X(:,1))
        SuffStats{2,i} = SuffStats{2,i} + all_X(j,:)*all_X(j,:)';
    end
    SuffStats{2,i} = SuffStats{2,i}/this.N_MC_samples;
    

    SuffStats{3,i} = mean(cell2mat(all_Y'),2);
    SuffStats{4,i} = 0;
    for j = 1: length(X(:,1))
        SuffStats{4,i} = SuffStats{4,i} + all_Y{j}*all_Y{j}';
    end
    SuffStats{4,i} =  SuffStats{4,i}/this.N_MC_samples;
    
%     SuffStats{5,i} = 0;
%     for j = 1: this.N_MC_samples
%         SuffStats{5,i} = SuffStats{5,i} + all_X(j,:)*all_X(j,:)'/this.N_MC_samples;
%     end
end


end










