classdef pc_lin_feature_gaussian1 < pc_lin_feature_gaussian
    %PCF_LIN_FEATURE_GAUSSIAN1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Init_var;
    end
    
    properties (Constant)
        N_param_coarse = 1;
    end
    
    methods
        % Constructor
        function this = pc_lin_feature_gaussian1(g,cg, c2f, feature_functions, init_var)
            
            this@pc_lin_feature_gaussian(g,cg, c2f, feature_functions);
            this.Init_var = init_var;
        end
        
        
        % Initialize optimal parameters
        function this = initParameters(this,varargin)
            if(numel(varargin) == 2)
                cov_matrix = varargin{1};
                coefficients = varargin{2};
            else
                % find out how many feature functions there are
                ncelem = this.Coarse_grid.N_elem;
                nfelem = this.Fine_grid.N_elem;
                testmat = this.Feature_functions(ones(1,nfelem));
                nfeat = length(testmat(1,:));
                
                cov_matrix = eye(this.N_param_coarse*ncelem)*this.Init_var;
                coefficients = ones(nfeat,1)/nfeat;
            end
            
            this.Cov_matrix = cov_matrix;
            this.Cov_matrix_chol = chol(cov_matrix);
            this.Det_cov_matrix = prod(diag(this.Cov_matrix_chol))^2;
            this.Log_det_cov_matrix = 2*sum(log(diag(this.Cov_matrix_chol)));
            this.Coefficients = coefficients;
        end
        
    end
    
end

