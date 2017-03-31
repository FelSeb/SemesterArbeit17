classdef EM
    %Expectation maximisation
    
    properties
        % Refinement PDF
        PCF;
        
        % Coarse graining PDF
        PC;
        
        % Primary model
        Primary_model;
        
        % Data directory (for now we store the training data directly)
        Data_provider;
            
        % Number of data points
        N_training_samples
        
        % Store sufficient statistics
        SuffStats;
        
        % Maximum number of iterations
        Max_iter;
        
        % Sammpler
        Sampler
    end
    
    methods
        
        % Constructor
        function this = EM(pc, pcf, primaryModel, n_training_samples, Data_provider, maxiter, sampler)
            this.PCF = pcf;
            this.PC = pc;
            this.Data_provider = Data_provider;
            
            this.Max_iter = maxiter;
            
            this.N_training_samples = n_training_samples;
            
            warning('primaryModel really needed?')
            this.Primary_model = primaryModel;
            
            this.Sampler = sampler;
        end
        
        % Return value of logQopt_i = log_pc_i + log_pcf_i
        [val,info] = logQopt_i(this,X,x_i,y_i, varargin)
        
        % Optimize parameters via EM scheme
        [theta_pc, theta_pcf, this] = optimizeEM(this);
        
        % Compute suffiecient statistics via MCMC
        SuffStats = getSuffStats(this);
        
        % Set training data
        function this = setData(x,y)
            this.x = x;
            this.y = y;
        end
        
        % This funciton returns training sample i. It loads training data 
        % batchwise into the workspace in order to economize storage space.
        function [x_i, y_i] = provideDataPoint(this,i)
            
            % Load index to batch map and its inverse
            if(isempty(this.I2B))
               load([this.Data_directory,'/I2B']);
               load([this.Data_directory,'/B2I']);
               this.I2B = I2B; 
               this.B2I = B2I; 
            end
            
            % load new batch
            if(isempty(this.current_batch) || i == 1 || this.I2B(1,i) > this.I2B(1,i-1))
                load([this.Data_directory,'/data_point_batch',num2str(this.I2B(1,i))])
                this.current_batch = data_point_batch;
            end
            % extract one data point from batch
            x_i = this.current_batch{this.I2B(2,i)}.x;
            y_i = this.current_batch{this.I2B(2,i)}.y;
        end
        
    end
    
%     methods (Static)
%         
%         function [ logp ] = gaussianLogProp( x,y, propvar )
%             % Gaussian proposal distribution. Probability density of choosing x given
%             % y. x and y are column vectors.
%             
%             k = length(x);
%             vec = (x-y);
%             logdet_cova = (1/propvar)^k;
%             logp = -k/2 * log(2*pi) - 0.5 * logdet_cova - 0.5 * (vec' *vec) / propvar;
%         end
%         
%         function [ smpl ] = gaussianPropRnd(y, propvar )
%             % Gaussian proposal distribution. Probability density of choosing x given
%             % y. x and y are column vectors.
%             smpl = mvnrnd(y',eye(length(y))*propvar)';
%         end
%         
%         function logp = MALALogProp(x, y, propvar, varargin)
%             if(numel(varargin) > 0)
%                info =  varargin{1};
%             else
%                info = cell(1,2);
%                info(:) = {0};
%                warning('MALALogProp not provided with gradient of logpdf')
%             end
%             
%             dt = 1;
%             sigma2 = propvar * dt;
%             mean = y' + sigma2/2 * info{1,2};
%             % cova = eye(length(y))*sigma2;
%             % logp = log(mvnpdf(x',mean,cova));
%             
%             k = length(mean);
%             vec = (x-mean');
%             logdet_cova = (1/sigma2)^k;
%             logp = -k/2 * log(2*pi) - 0.5 * logdet_cova - 0.5 * (vec' *vec) / sigma2;
%         end
%         
%         function smpl = MALAPropRnd(y, propvar, varargin)
%             if(numel(varargin) > 0)
%                info =  varargin{1};
%             else
%                info = cell(1,2);
%                info(:) = {0};
%                warning('MALALogProp not provided with gradient of logpdf')
%             end
%             
%             dt = 1;
%             sigma2 = propvar * dt;
%             mean = y' + sigma2/2 * info{1,2};
%             cova = eye(length(y))*sigma2;
%             
%             smpl = mvnrnd(mean,cova)';
%         end
%         
%         
%         % Modified matlab Metropolis-Hastings sampler
%         [smpl,accept] = MYmhsample(start,nsamples,varargin)
%         
%     end
    
end

