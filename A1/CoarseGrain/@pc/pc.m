classdef pc
    %PC Summary of this class goes here
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
        
        % Fine grained material parameter data. This is a struct with as
        % many entries as there are training samples. In a later
        % implementation this may be replaced with a directory containing 
        % batches of training data. This becomes necessary if there are
        % many training samples of high dimension.
        x;
        
        
    end
    
    methods
        % Constructor
        function this = pc(fine_grid, coarse_grid)
            this.Fine_grid = fine_grid;
            this.Coarse_grid = coarse_grid;
        end
    end
    
    methods
        % Setting vector of FG values
        function this = setFineVals(this,y)
            this.y =  y;
        end
        
    end
    
    methods (Abstract)
%         % Coarse to fine mapping
%         coarse_field = fine2coarse(this, fine_field);
        
        % Evaluate PDF
        pdf_val = evalPDF(this, X, x);
        
        % Evaluate log PDF
        [log_pdf_val, grad_log_pdf_val] = evalLogPDF(this, X, x, varargin);
        
    end
    
end

