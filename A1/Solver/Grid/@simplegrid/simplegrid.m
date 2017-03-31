classdef simplegrid
    %GRID Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Number of elements
        N_elem;
        
        % Number of nodes
        N_node;
        
        Node2eq; % Map: node -> equation
        Elem2node; % Map: element -> node
        Node2coord; % Map: node -> coordinate
        
        % Node numbers of all non-dirichlet-boundary nodes
        All_nondir_nodes;
        % Node numbers of all boundary nodes
        All_boundary_nodes;
        % Element numbers of all boundary elems
        All_boundary_elems;
        % Node numbers of Dirichlet boundary nodes (one array for every edge)
        Dir_boundary_nodes;
        % Node numbers of Dirichlet boundary elements (one array for all)
        Dir_boundary_elems;
        % Node numbers of v. Neumann boundary nodes (one array for all)
        Neu_boundary_nodes;
        % Node numbers of Dirichlet boundary nodes (one array for every edge)
        Neu_boundary_elems;
        
        % Material parameters are assumed to be constant troughout an
        % element. The corresponding map (element number -> parameter) is
        % Elem2param. It is a matrix of dimension 2 x 2 x N_elem:
        Elem2param;
        
        % Material parameters in vector format. It takes one of the three
        % following forms (Xelemnumber_paramnumber):
        % [X1_1, X1_2, X1_3, X2_1, X2_2,... ]'
        % [X1_1, X1_2, X2_1, X2_2, X3_1,... ]'
        % [X1_1, X2_1, X3_1, X4_1, X5_1,... ]'
        % X is parametrized as follows (ensureing posiive definiteness of 
        % lambda):
        % lambda = L'*L;  L = [exp(X_1), X_3; 0, exp(X_2)];
        X;
        
        % Finiete element type
        FE;
    end
    
    methods
        % set Elem2param map from parameter vector X
        function this = setElem2param(this,X, paramTrafo)
            E2Pdim = size(this.Elem2param);
            if(length(E2Pdim) == 3)
                nparam = 3;
            elseif(E2Pdim(1) == 2)
                nparam = 2;
            else
                nparam = 1;
            end
            this.Elem2param = paramTrafo.getElem2ParamfromX(X,nparam);
        end
        
        % set unconverted parameter vector X 
        function this = setX(this, X)
            this.X = X;
        end
        
        % convert Elem2param map to parameter voector X
        function this = makeXfromElem2Param(this, paramTrafo)
            this.X = paramTrafo.getXfromElem2param(this.Elem2param);
        end
        
    end
    
    methods (Abstract)
        % Method to derive maps from the grid strucure: Elem2node,
        % Node2coord, All_boundary_elems, All_boundary_nodes, Elem2boundary
        this = makeMaps(this);
        
        
    end
    
end

