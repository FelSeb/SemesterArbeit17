classdef pcf
    %PCF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % x and y nodal coordinates of the high resolution finite element
        % mesh. Format: Matricies with as many rows as nodes per element
        % and as many columns as elements (nodes_per_elem x N_elem x dim).
        Fine_grid;
        
        % Map telling which fine grid elements belong to which coarse grid
        % element. Format: Each column number corresponds to a coarse
        % element number. The rows contain the numbers of the corresponding
        % fine elements
        C2F;
        
        % x and y nodal coordinates of the low resolution finite element
        % mesh. Format: Matricies with as many rows as nodes per element
        % and as many columns as elements (nodes_per_elem x N_elem x dim).
        Coarse_grid;
        
        % Solver object that allows to solve the PDE on the coarse grid
        Solver
        
        % Number of coarse and fine elements
        N_coarse;
        N_fine;
    end
    
    methods
        % Constructor
        function this = pcf(fine_grid, coarse_grid, c2f, solver)
            this.Fine_grid = fine_grid;
            this.Coarse_grid = coarse_grid;
            this.C2F = c2f;
            this.Solver = solver;
            % Set_up solver:
            % 1.) Set grid
            this.Solver = this.Solver.setGrid(this.Coarse_grid);
            % 2.) Preevaluate shape- /test- functions (and derivatives)
            this.Solver = this.Solver.preEval();
        end
    end
    
    methods (Abstract)
        % Coarse to fine mapping
        fine_field = coarse2fine(this, coarse_field);
        
        % Evaluate PDF
        [pdf_val, Y] = evalPDF(this, y, Y);
        
        % Evaluate log PDF
        [log_pdf_val, Y, grad_log_pdf_val] = evalLogPDF(this, y, Y, varargin);
        
    end
    
end

