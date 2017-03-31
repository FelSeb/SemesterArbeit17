function [ SuffStats ] = getSuffStats( this )
% Computing a MCMC Method to compute the sufficient statistics.

n_coarse_elem = this.PCF.Coarse_grid.N_elem;

E2Pdim = size(this.PCF.Coarse_grid.Elem2param);

% Proposal distribution. Probability density of choosing x given y

LOGPROPPDF = @(X_new,X_old, varargin) this.Proposer.logproppdf( X_new, X_old, varargin{:});

% Random number generator of proposal distribution:
PROPRND = @(X_old, varargin) this.Proposer.proprnd(X_old, varargin{:});

% This for loop may be parallelized!!!
for i = 1:this.N_training_samples
    % set training data x and y of current batch
    [x_i, y_i] = this.Data_provider.provideDataPoint(i);
    
    % Target distribuiton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LOGPDF = @(X) logQopt_i(this,X,x_i,y_i,  1  );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set initial guess. Using the mean of pc seems to be a convenient
    % choice.
    init = (this.PC.Feature_functions(x_i) * this.PC.Coefficients)';
    
    % Determine appropriate stepwidth

    tic
    % Generate all X-samples
    [smpl, accept] = MYmhsample(init,this.N_MC_samples,...
        'logpdf',LOGPDF,...
        'logproppdf',LOGPROPPDF,'proprnd',PROPRND,...
        'returninfo',1);
    toc
    
    % Derive <X>, <|X|^2>, <Y(X)>, <Y(X)Y(X)'>
    %all_X = smpl;
    all_X = smpl.samples;
    all_Y = smpl.info;
    
    %%%%%%%%%%%%%%%%%%% for test purposes &%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dvals = smpl.info2;
    auto_cov = estimateAutoCovariance(all_X,1/5);
    accept
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SuffStats{1,i} = mean(all_X',2);
    SuffStats{2,i} = 0;
    for j = 1: this.N_MC_samples
        SuffStats{2,i} = SuffStats{2,i} + all_X(j,:)*all_X(j,:)';
    end
    SuffStats{2,i} = SuffStats{2,i}/this.N_MC_samples;
    

    SuffStats{3,i} = mean(cell2mat(all_Y'),2);
    SuffStats{4,i} = 0;
    for j = 1: this.N_MC_samples
        SuffStats{4,i} = SuffStats{4,i} + all_Y{j}*all_Y{j}';
    end
    SuffStats{4,i} =  SuffStats{4,i}/this.N_MC_samples;
    
%     SuffStats{5,i} = 0;
%     for j = 1: this.N_MC_samples
%         SuffStats{5,i} = SuffStats{5,i} + all_X(j,:)*all_X(j,:)'/this.N_MC_samples;
%     end
end


end

% 
% function rho = estimateAutoCovariance(all_X,jfrac)
% [nsamples,ndim] = size(all_X);
% 
% % sample mean
% mu_hat = mean(all_X,1);
% 
% % sample covariance
% Sigma_hat = 0;
% for n = 1:nsamples
%     Sigma_hat = Sigma_hat + (all_X(n,:)-mu_hat).^2;
% end
% Sigma_hat = Sigma_hat/nsamples;
% 
% nj = ceil(nsamples*jfrac);
% rho = zeros(nj+1,ndim);
% for j = 0:nj
%     for n = 1:nsamples-j
%         rho(j+1,:) = rho(j+1,:) + (all_X(n,:)-mu_hat).*(all_X(n+j,:)-mu_hat);
%     end
%     rho(j+1,:) = 1/(nsamples-j)*rho(j+1,:)./Sigma_hat;
% end
% figure;
% plot(rho);
% % figure;
% % plot(all_X);
% % figure
% % plot(cummean(all_X))
% 
% end


function stepWidth = optimizeStepwidth(hi, lo, opt, init, LOGPDF, LOGPROPPDF, PROPRND)

lo = 0.5;
hi = 2;
opt = 0.4;

c = lo;
b = ((1-lo) -(hi-lo)*opt^2)/(opt-opt^2);
a = hi-lo-b;




fac = a*x.^2 + b*x + c;

N_test_samples = 111;


% Generate all X-samples
[smpl, accept] = MYmhsample(init,N_test_samples,...
    'logpdf',LOGPDF,...
    'logproppdf',LOGPROPPDF,'proprnd',PROPRND,...
    'returninfo',1)

fac = a*x.^2 + b*x + c;

end






