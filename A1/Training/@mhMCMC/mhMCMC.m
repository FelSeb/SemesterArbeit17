classdef mhMCMC
    %MCMC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Proposer
        
        Logtarget
        
        N_MC_samples
        
        LOGPROPPDF
        
        PROPRND
    end
    
    methods
        % Constructor
        function this = mhMCMC(proposer, n_mc_samples)
            if(~isa(proposer,'PropDensity'))
                error('Proposer must be an object of type PropDensity')
            end
            
            this.Proposer = proposer;
            
            this.N_MC_samples = n_mc_samples;
            
            this.LOGPROPPDF = @(X_new,X_old, varargin) this.Proposer.logproppdf( X_new, X_old, varargin{:});
            this.PROPRND = @(X_old, varargin) this.Proposer.proprnd(X_old, varargin{:});
            
        end
        
        
        function this = setProposalSteplength(prop_var)
            this.Proposer.Prop_step_length = prop_var;
        end
        
        function [samples, accept] = generateSamples(this,LOGPDF, init)
            [samples,accept] = this.MYmhsample(init,this.N_MC_samples,...
                'logpdf',LOGPDF,...
                'logproppdf',this.LOGPROPPDF,'proprnd',this.PROPRND,...
                'returninfo',1);
        end
        
        
        function this = optimizeStepwidth(this,init, LOGPDF)
            nmax_attemt = 5;
            N_test_samples = 111;
            lo_limit = 0.3;
            hi_limit = 0.8;
            lo = 0.5;
            hi = 2;
            opt = 0.4;
            
            c = lo;
            b = ((1-lo) -(hi-lo)*opt^2)/(opt-opt^2);
            a = hi-lo-b;
            
            n_attempt = 1;
            accept = -1;
            while((accept < lo_limit || accept > hi_limit) && n_attempt < nmax_attemt)
                % Generate all X-samples
                [smpl, accept] = this.MYmhsample(init,N_test_samples,...
                    'logpdf',LOGPDF,...
                    'logproppdf',this.LOGPROPPDF,'proprnd',this.PROPRND,...
                    'returninfo',1);
                
                fac = a*accept^2 + b*accept + c;
                this.Proposer.Prop_step_length = this.Proposer.Prop_step_length * fac;
            end
            
        end
        
    end
    
    
    methods(Static)
        % Metropolis hastings algorithm
        [smpl,accept] = MYmhsample(start, nsamples, varargin)
        
        
        function rho = estimateAutoCovariance(all_X,jfrac)
            [nsamples,ndim] = size(all_X);
            
            % sample mean
            mu_hat = mean(all_X,1);
            
            % sample covariance
            Sigma_hat = 0;
            for n = 1:nsamples
                Sigma_hat = Sigma_hat + (all_X(n,:)-mu_hat).^2;
            end
            Sigma_hat = Sigma_hat/nsamples;
            
            nj = ceil(nsamples*jfrac);
            rho = zeros(nj+1,ndim);
            for j = 0:nj
                for n = 1:nsamples-j
                    rho(j+1,:) = rho(j+1,:) + (all_X(n,:)-mu_hat).*(all_X(n+j,:)-mu_hat);
                end
                rho(j+1,:) = 1/(nsamples-j)*rho(j+1,:)./Sigma_hat;
            end
            figure;
            plot(rho);
            % figure;
            % plot(all_X);
            % figure
            % plot(cummean(all_X))
            
        end
    end
    
end

