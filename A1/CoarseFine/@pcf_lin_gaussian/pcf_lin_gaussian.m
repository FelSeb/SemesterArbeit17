classdef pcf_lin_gaussian < pcf
    
    % y = mu + W * Y + z; z ~ N(0,S)
    
    properties
        % Linear transformation matrix for fast mapping of nodal values of
        % a single coarse to the nodal values of the corresponding fine
        % elements. It is of dimension (dim_fine x dim_coarse). If not all
        % coarse elements have the same shape there must be as many
        % interpolation matricies as there are elements.
        Interpolation_matrix;
        
        % Offset-vector
        Offset;
        
        % Covariance matrix, its Cholesky factorization and determinant
        S;
        S_chol;
        Det_S;
        Log_det_S;
        
        % String specifying which parameters are to be optimized and which
        % are not.
        Var_param_string;
        
        % Boolean specifying whether all elements have the same shape
        Isoshape;
        
        % initial variance
        Init_var
        
    end
    
    properties (Constant)
        % Is there a function for a non-iterative computation of an optimum
        % of the lower bound (see EM)?
        Direct_optim = 1;
    end
    
    methods
        % Constructor
        function this = pcf_lin_gaussian(fg,cg,c2f,variable_params_string, solver, isoshape, init_var)
            % Constructor of super class
            this@pcf(fg, cg, c2f, solver);
            
            % in case no optimisation of the parameters is to be performed
            this.Var_param_string = variable_params_string;
            
            warning('is isoshape really needed?')
            this.Isoshape = isoshape;
            
            this.Init_var = init_var;
        end
        
        % Compute Intetpolation matrix
        this = getInterpMat(this, isoshape);
        
        
        % Bilinear interpolation for refinement:
        fine_field = bilinQuadRefine(this, coarse_field )
        
        % Overriede method:
        function fine_field = coarse2fine(this, coarse_field)
            fine_field = this.bilinQuadRefine(coarse_field);
        end
        
        
        % Evaluate PDF or log PDF
        [pdf_val, Y, grad_pdf_val] = evalGaussianPDF(this, y, X, eval_log, varargin)
        
        % Override Method: Evaluate PDF
        function [pdf_val,y] = evalPDF(this, y, X)
            pdf_val = this.evalGaussianPDF(y,X,0);
        end
        
        % Override Method: Evaluate log PDF
        function [log_pdf_val, y, grad_log_pdf_val] = evalLogPDF(this, y, X, varargin)
            if(numel(varargin)>0)
                compute_grad = varargin{1};
            else
                compute_grad = 0;
            end
            
            if(~compute_grad)
                [log_pdf_val, y] = this.evalGaussianPDF(y, X, 1);
                grad_log_pdf_val = [];
            else
                [log_pdf_val, y, grad_log_pdf_val] = this.evalGaussianPDF(y, X, 1, 1);
            end
        end
        
        % Compute optimal parameters given the sufficient statistics
        function this = getOptimParams(this, suff_stats, data_provider, n_training_samples)
            if(strcmp(this.Var_param_string,'allfix'))
                % Do nothing;
            elseif(strcmp(this.Var_param_string,'varS'))
                this = this.getOptimSDiag(suff_stats, data_provider, n_training_samples);
            else
                error('Optimisation routine not implemented yet')
            end
        end
        
        % Compute optimal noise covariance matrix
        this = getOptimCovMatrix(this);
        
        
        % Initialize optimal parameters
        function this = initParameters(this,varargin)
            if(numel(varargin) == 2)
                this.S = varargin{1};
                this.Offset = varargin{2};
            else
                this.S = eye(this.Fine_grid.N_node) * this.Init_var;
                this.Offset = zeros(this.Fine_grid.N_node,1);
            end
            
            this = this.getInterpMat;
            
            this.S_chol = chol(this.S);
            this.Det_S = prod(diag(this.S_chol))^2;
            this.Log_det_S = 2*sum(log(diag(this.S_chol)));
        end
    end
    
end

