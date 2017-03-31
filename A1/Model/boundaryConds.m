classdef boundaryConds
    %BOUNDARY_CONDS 
    
    properties
        % A global function from which boundary conditions are derived
        Global_fun;
        
        % The domain
        Domain;
        
        % This is a cell array. There must be as many cell columns as there
        % are edges. And two function-handels per column (Dir. and Neu.)
        Bound_cond_funs;
        
        Info;
    end
    
    methods
        % Constructor
        function this = boundaryConds(domain)
            this.Domain = domain;
        end
        
        % Evaluate boundary functions of an edge
        function this = SetUpPolynimialBoundConds(this, boundaryCoeffs)
            
            this.Bound_cond_funs = cell(2,this.Domain.N_edges);
            for bc = 1:this.Domain.N_edges
                % Temperature field generating Dirichlet boundary conditions
                this.Bound_cond_funs{1,bc} = @(x,y)...
                    boundaryCoeffs(1) + boundaryCoeffs(2)*x + boundaryCoeffs(3)*y + boundaryCoeffs(4)*x.*y;
                
                % Temperature gradient generating v. Neumann boundary conditions
                this.Bound_cond_funs{2,bc} = @(x,y)...
                    [boundaryCoeffs(2)+ boundaryCoeffs(4)*y; boundaryCoeffs(3) + boundaryCoeffs(4)*x];
            end
            
            this.Info = struct;
            this.Info.Coeffs = boundaryCoeffs;
        end
        
        
        function reference_sol = evaluateBoundaryFun(this, x_grid, y_grid, bc)
            reference_sol = this.Bound_cond_funs{1,bc}(x_grid,y_grid);
        end
    end
    
end

