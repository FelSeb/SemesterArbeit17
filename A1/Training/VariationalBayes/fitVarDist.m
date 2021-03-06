function [ varDist ] = fitVarDist(ptarget, varDist)
% ptarget is a function handle representing the (unnnormalized target pdf).
% It must be of the following form [val,info] = ptarget(X), where info is a
% cell array and costructed as follows:
%     info{1,1} = Y(X); (Coarse solution)
%     info{1,2} = dval; (gradient of the target distribution w.r.t X)
% val is the value of the target distribution at input X.


N_dim = varDist.N_dim;

standard_mean = zeros(1,N_dim);
standard_cov = eye(N_dim);


% Good settings for diag gaussian: N_maxiter = ca.1000; N_maxmc_sample = 1,
% eta = 0.1;
% Good settings for full gaussian: N_maxiter = ca.500; N_maxmc_sample = 10,
% eta = 0.1; (Not better in terms of ELBO value)


N_maxiter = 800;
N_maxmc_sample = 1;

% Adam parameters:
eta = 0.1;
beta1 = 0.9;
beta2 = 0.999;
smalleps = 1e-8;
m = 0;
v = 0;

%ELBO = zeros(1,N_maxiter);
%mu = zeros(N_maxiter, N_dim);
%var_diag = zeros(N_maxiter, N_dim);

% Maximization of the eveidence lower bound:
for iter = 1:N_maxiter

    % Compute gradient of entropy term @ phi_current w.r.t. phi
    DEntropyDphi = varDist.getGradEntropy_dphi;
    %Entropy = varDist.getEntropy;
    
    % sample N_maxmc_sample epsilon from the standard Normal distribution
    epsilon = mvnrnd(standard_mean, standard_cov, N_maxmc_sample);
    
    % Initialize MC average of the noisy gradient
    mc_av_grad = zeros(1, varDist.N_var_param);
    %mc_av_val = 0;
    
    % Generate MC samples
    for mc_sample = 1:N_maxmc_sample
        % Transform epsilon to X
        X = varDist.reparamFun(epsilon(mc_sample,:));
        
        % Evaluate grad_X(ptarget(X( phi_current, epsilon(mc_sample) )))
        [val,info] = ptarget(X);
        DTargetDX = info{1,2};
        
        % Evaluate grad_phi(X( phi_current,  epsilon(mc_sample) ))
        dXdphi = varDist.gradReparamFun_dphi(epsilon(mc_sample,:));
        
        % Update MC average of the noisy gradient
        mc_av_grad  = mc_av_grad + DTargetDX * dXdphi;
        %mc_av_val = mc_av_val + val;
    end
    mc_av_grad = mc_av_grad/N_maxmc_sample;
    %mc_av_val = mc_av_val/N_maxmc_sample;
    
    %elbo = mc_av_val + Entropy;
    %ELBO(iter) = elbo;
    
    % Add entropy gradient and Expectation gradient
    complete_grad = DEntropyDphi + mc_av_grad;
    
    % Apply one of the stochastic gradient update rules to find new
    % phi_current. (AdaGrad or Adam or ...)
    
    % Adam (Adaptive moment estimation)
    m = beta1*m + (1-beta1)*complete_grad;
    v = beta2*v + (1-beta2)*complete_grad.^2;
    
    m_hat = m/(1-beta1^iter);
    v_hat = v/(1-beta2^iter);
    
    delta_phi = eta * m_hat./(sqrt(v_hat) + smalleps);
    
    varDist = varDist.updatePhi(delta_phi);
    
    %mu(iter,:) = varDist.Phi_mu;
    %var_diag(iter,:) = diag(varDist.Sigma_chol' * varDist.Sigma_chol);
end


    
% plot(ELBO);
% title('ELBO')
% figure;
% plot(mu)
% title('MU')
% figure;
% plot(var_diag);
% title('VARDIAG')

% Return fitted variational distribution
end








