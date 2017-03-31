classdef polyangularDomain
    
    
    properties
        Domain_ID;
        
        % Defines the edges of the domain. Each colum conatains the
        % coordinates of the two nodes of an edge. Neighbouring columns
        % corrspond to edges that have a common node.
        Edges;
        
        % Define the smalles rectangle enclosing the domain
        X_min;
        X_max;
        Y_min;
        Y_max;
        
        % Number of edges
        N_edges;
        
        % Normal vectors for all edges defined in domain
        All_edges_normvec;
    end
    
    methods
        
        % Compute vectors that are normal to all edges of the domain
        function this = GetNormVecs(this)
            % Get the outward pointing normal vectors for each edge of the domain
            this.All_edges_normvec = zeros(2,this.N_edges);
            for e = 1:this.N_edges
                vec = this.Edges([3,4],e) - this.Edges([1,2],e);
                
                this.All_edges_normvec(1,e) = vec(2);
                this.All_edges_normvec(2,e) = -vec(1);
                
                this.All_edges_normvec(:,e) = this.All_edges_normvec(:,e)/norm(this.All_edges_normvec(:,e));
            end
        end
        
    end
end

