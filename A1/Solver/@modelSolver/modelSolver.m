classdef modelSolver
    %MODELSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % primaryModel object
        PrimaryModel;
        
        Grid; % Discretized domain
        
        % Preevaluated shape- /test- functions (and derivatives) at Gauss
        % points
        N_all;
        DN_all;
        N_all1D;
        DN_all1D;
        
        % Stiffness matrix C
        C;
        
        % Gradient of the stiffness matrix
        C_grad;
        
        % Rhight hand side
        Rhs;
        
        % Gradient of the Rhs vector
        Rhs_grad;
        
        % Solution vector
        Sol;
        
        % Gradient of the solution vector w.r.t. material parameters
        % Bad idea:
        Sol_grad
        
        % Postprocessed solution
        PostSol;
        
        % Transformation from lambda to X (and inverse)
        Param_trafo;
    end
    
    methods
        % Constructor
        function this = modelSolver(a_primaryModel, varargin)
            this.PrimaryModel = a_primaryModel;
            if(numel(varargin) > 0 && isa(varargin{1}, 'paramTrafo'))
                this.Param_trafo = varargin{1};
            elseif(numel(varargin) > 0 )
               error('Object of type paramTrafo required as second input argument') 
            end
        end
        
        % Preevaluate shape- /test- functions (and derivatives)
        this = preEval(this);
        
        % Compose stiffness matrix
        this = composeStiffness(this, varargin);
        
        % Apply Boundary conditions
        this = applyBoundCond(this, varargin);
        
%         % Apply Boundary conditions
%         this = applyBoundCondOLD2(this, varargin);
        
        % Solve the system
        this = solveSystem(this, varargin);
        
        % Postprocess solution
        this = getPostSol(this, Compute_gradient, varargin);
        
        % Set grid
        function this = setGrid(this, grid)
            this.Grid = grid;
        end
        
        % Plot the solution
        function plotSol(this, Display_gradient)
            this.Grid.plotSol(this.PostSol, Display_gradient)
        end
        
        % Method that returns dirichlet boundary values and other things in
        % a convenient format [dirichlet equation numbers, corresponding
        % dirichlet values, non-dirichlet equation numbers]
        [dirs, dir_vals, inner] = formatDirBoundInfo(this);
        
        % Evaluate model
        function [res, varargout] = evaluateModel(this, varargin)
            %Usage:
            % [Y, this, dYdX] = evaluateModel(this, X, computeGrad, computeSolGrad)
            % The last three input arguments are optional. If the model is to
            % be evaluated for a given parameter field, insert the
            % corresponding field in vector format at the placeholder X;
            % else: insert any scalar. If additionally the gradient of the
            % FEM matrix and right hand side w.r.t. the parameter field is 
            % to be computed set computeSolGrad to one; else: zero.
            % If additionally the gradient of the solution vector w.r.t. 
            % the parameter field is to be computed 
            % set computeSolGrad to one; else: zero.
            
            if(numel(varargin) > 1 && numel(varargin) < 4)
                % Compute gradient of solution vector w.r.t. parameter
                % field if varargin{2} == 1.
                if(numel(varargin{1}) > 0)
                    % Compute output for a given parameter field.
                    this.Grid = this.Grid.setElem2param(varargin{1},this.Param_trafo);
                    this.Grid = this.Grid.setX(varargin{1});
                    this = this.composeStiffness(varargin{2});
                else
                    % compute output for a randomly drawn parameter field.
                    this.Grid = this.Grid.remakeElem2Param(this.PrimaryModel);
                    if(varargin{2})
                        this.Grid = this.Grid.makeXfromElem2Param(this.Param_trafo);
                        warning('Is this really necessary?') 
                    end
                    this = this.composeStiffness(varargin{2});
                end
                this = this.applyBoundCond(varargin{2});
                
                if(numel(varargin) == 3 && varargin{3} && varargin{2})
                    this = this.solveSystem(varargin{3});
                    this = this.getPostSol(0,varargin{3});
                elseif(numel(varargin) == 3 && varargin{3} && ~varargin{2})
                    error('If computeSolGrad == 1, computeGrad hast to be equal to 1 aswell!')
                else
                    this = this.solveSystem();
                    this = this.getPostSol(0);
                end
            elseif(numel(varargin) < 2)
                if(numel(varargin) == 1 && numel(varargin{1}) > 1)
                    this.Grid = this.Grid.setElem2param(varargin{1}, this.Param_trafo);
                    this = this.composeStiffness();
                else
                    this.Grid = this.Grid.remakeElem2Param(this.PrimaryModel);
                    this = this.composeStiffness();
                end
                this = this.applyBoundCond();
                this = this.solveSystem();
                this = this.getPostSol(0);     
            else
                error('Number of input arguments not admissible')
            end
            
            res = this.PostSol.T;
            
            if(nargout > 1)
                varargout{1} = this;
            end
            if(nargout > 2)
                varargout{2} = this.PostSol.dYdX;
            end
            
        end
        
    end
    
    
end

