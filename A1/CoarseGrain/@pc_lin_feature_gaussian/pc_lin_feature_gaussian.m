classdef pc_lin_feature_gaussian < pc
    
    % y = mu + W * Y + z; z ~ N(0,S)
    
    properties
        % Coefficients of the linear combination of feature functions
        Coefficients;
        
        % Function handle of function returning values of all feature
        % function in a 2D-array (each column one feature).
        Feature_functions;
        
        % Covariance matrix, its Cholesky factorization and determinant
        Cov_matrix;
        Cov_matrix_chol;
        Det_cov_matrix;
        Log_det_cov_matrix;
        
        % Number of coarse elements
        N_elem_coarse;
        
        % Coarse to fine mapping
        C2f;

    end
    
    properties (Constant)
        % Is there a function for a non-iterative computation of an optimum
        % of the lower bound (see EM)?
        Direct_optim = 1;
    end
    
    methods
        % Constructor
        function this = pc_lin_feature_gaussian(fg,cg, c2f, feature_functions)
            % Constructor of super class
            this@pc(fg, cg);
            this.Feature_functions = feature_functions;
            this.C2f = c2f;
        end
        
        %         % Overriede method:
        %         function coarse_field = fine2coarse(this, fine_field)
        %             coarse_field = this.featureCoarseGrain(fine_field);
        %         end
        %
        %         % Coarsegraining funciton
        %         coarse_field = featureCoarseGrain(this, fine_field)
        
        % Evaluate PDF or log PDF
        [pdf_val, grad_pdf_val] = evalGaussianPDF(this, X, x, eval_log, varargin)
        
        % Override Method: Evaluate PDF
        function pdf_val = evalPDF(this, X, x)
            pdf_val = this.evalGaussianPDF(X, x,0);
        end
        
        % Override Method: Evaluate log PDF
        function [log_pdf_val, grad_log_pdf_val] = evalLogPDF(this, X, x, varargin)
            if(numel(varargin)>0)
                compute_grad = varargin{1};
            else
                compute_grad = 0;
            end
            
            if(~compute_grad)
                [log_pdf_val, grad_log_pdf_val] = this.evalGaussianPDF(X, x, 1, compute_grad); 
            else
                [log_pdf_val, grad_log_pdf_val] = this.evalGaussianPDF(X, x, 1, compute_grad);
            end
            
        end
        
        % Compute optimal parameters given the sufficient statistics and no
        % prior
        this = getOptimParamsNoPrior(this, suffStats, data_provider, n_training_samples)
        
        % Compute optimal parameters given the sufficient statistics
        function this = getOptimParams(this, suffStats, data_provider, n_training_samples)
            %this = getOptimParamsNoPrior(this, suffStats, data_provider, n_training_samples);
            this = getOptimParamsSparsity(this, suffStats, data_provider, n_training_samples);
        end
        
    end
    
    methods (Abstract)
        % Initialize optimal parameters
        this = initParameters(this,varargin)
    end
    
end

