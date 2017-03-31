function [ SuffStats ] = getSuffStatsVarBay( this )
% Computing a MCMC Method to compute the sufficient statistics.

% This for loop may be parallelized!!!
SuffStats = cell(4,this.N_training_samples);
tic;
parfor i = 1:this.N_training_samples % (in parallel)
    % set training data x and y of current batch
    [x_i, y_i] = this.Data_provider.provideDataPoint(i);
    
    % Target distribuiton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LOGPDF = @(X) logQopt_i(this, X, x_i, y_i, 1);
    LOGPDF = @(X) this.logQopt_i(X, x_i, y_i, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Mean of pc:
    pcmean = (this.PC.Feature_functions(x_i) * this.PC.Coefficients)';
    phi_mu_init = pcmean;%log([1,1,1,1])/2;
    
    
    % In the case of varDistGaussian the format of phi_L_init determines
    % which type of gaussian is used as a variational distribution.
    % (Isotropic, orthotropic or full covariance matrix)
    %
    % orthotropic covariance matrix:
    N_dim = length(pcmean);
    % The initial variance should be chosen carefully: Let's assume 1 <
    % Lambda < 100. Then, since lambda = exp(X * 2) ->  0 < X < 2.3 ->
    % standard deviation is roughly 0.5 -> sigma2 = 0.25
    sigma2_init = 0.25;
    Sigma_diag_init = ones(1,N_dim)*sigma2_init;
    phi_L_init = log(Sigma_diag_init)/2;
    % full covariance matrix:
    % phi_L_init = [phi_L_init,ones(1,(N_dim-1)*N_dim/2)*0];
    
    varDist = varDistGaussian(phi_mu_init, phi_L_init);
    
    
    varDist = fitVarDist(LOGPDF, varDist);
    
    % Generate X samples
    N_MC_samples = 1000;
    X_samples = varDist.generateSamples(N_MC_samples);
    
    % Generate Y samples
    Y_samples = zeros(N_MC_samples, this.PC.Coarse_grid.N_node);
    for si = 1:N_MC_samples
        Y_samples(si, :) = this.PCF.Solver.evaluateModel(X_samples(si,:));
    end
    
    %%%%%%%%%% Visualtization of target and fitted distribution%%%%%%%%%%%%
    %     Sigma = varDist.Sigma_chol * varDist.Sigma_chol;
    %     F_gaussian = @(X) testTarget(X',varDist.Phi_mu,Sigma);
    %     plot1DDistribution(F_gaussian,varDist.Phi_mu);
    %     title('approx')
    %
    %     F_target = @(X) LOGPDF(X);
    %     plot1DDistribution(F_target,varDist.Phi_mu);
    %     title('ref')
    %%%%%%%%%%%%%%%%%%%%%% Test target distribuiton%%%%%%%%%%%%%%%%%%%%%%%%
    %     LOGPDF = @(X) testTarget(X, [0,0,0,0], eye(4));
    %     phi_mu_init = [3,3,3,3];
    %     Sigma_diag_init = [3,3,3,3];
    %     phi_L_init = log(Sigma_diag_init)/2;
    %
    %     varDist = varDistGaussian(phi_mu_init, phi_L_init);
    %
    %     varDist = fitVarDist(LOGPDF, varDist);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    susta1 = varDist.Phi_mu;
    susta2 = 0;
    for j = 1: length(X_samples(:,1))
        susta2 = susta2 + X_samples(j,:)*X_samples(j,:)';
    end
    susta2 = susta2/N_MC_samples;    

    susta3 = mean(Y_samples,1);
    susta4 = 0;
    for j = 1: length(X_samples(:,1))
        susta4 = susta4 + Y_samples(j,:)'*Y_samples(j,:);
    end
    susta4 =  susta4/N_MC_samples;

    SuffStats(:,i) = {susta1; susta2; susta3; susta4};   
end
toc;



end



function plot1DDistribution(F,Xloc)
n = 200;
X = linspace(-2,3,n);
Xfix = repmat(Xloc(1:end-1)',1,n);
X = [Xfix; X];

logp=zeros(n,1);

for i = 1:n
    logp(i) = F(X(:,i));
end

normalization = -max(logp);

p = exp(logp+normalization);

figure;
plot(X(end,:),p,'x')

end






