function this = makeMaps(this, PrimaryModel)
            % Map: element -> node
            this.Elem2node = this.FE.GetElem2node(this);
            
            % Node numbers of all boundary nodes
            this.All_boundary_nodes = [this.Topo_node(1:end-1,1)',...   % a
                this.Topo_node(end,1:end-1),...                         % b
                this.Topo_node(2:end,end)',...                          % c
                this.Topo_node(1,2:end)];                               % d
            
            % Element numbers of all boundary elems
            this.All_boundary_elems = [this.Topo_elem(1:end-1,1)',...   % a
                this.Topo_elem(end,1:end-1),...                         % b
                this.Topo_elem(2:end,end)',...                          % c
                this.Topo_elem(1,2:end)];                               % d
            
            % Node numbers of Dirichlet boundary nodes
            this.Dir_boundary_nodes{1} = [1];
            this.Dir_boundary_nodes{2} = [];
            this.Dir_boundary_nodes{3} = [];
            this.Dir_boundary_nodes{4} = [];
            
            All_Dir_boundary_nodes = [];
            for edge = 1:numel(this.Dir_boundary_nodes)
                All_Dir_boundary_nodes = [All_Dir_boundary_nodes,this.Dir_boundary_nodes{edge}];
            end
            
            % Node numbers of v. Neumann boundary nodes. (Not needed?)
            this.Neu_boundary_nodes = setdiff(this.All_boundary_nodes, All_Dir_boundary_nodes);
            
            
            % Cellarray of dimansion (Number of edges).
            % {[elem_1_of_edge_1, elem_2_of_edge_1,...], [elem_1_of_edge_2, elem_2_of_edge_2,...],...}
            % Number of cell entries = Number of edges
            nbe{1} = this.Topo_elem(:,1)';   % a
            nbe{2} = this.Topo_elem(end,:);  % b
            nbe{3} = this.Topo_elem(:,end)'; % c
            nbe{4} = this.Topo_elem(1,:);   % d
            
            % Map that indicates which nodes of an element belong to the
            % boundary
            nbes{1} = repmat([1;2],[1,numel(nbe{1})]);
            nbes{2} = repmat([2;3],[1,numel(nbe{2})]);
            nbes{3} = repmat([4;3],[1,numel(nbe{3})]);
            nbes{4} = repmat([1;4],[1,numel(nbe{4})]);
            
            % Compose into a single class attribute
            this.Neu_boundary_elems = {nbe,nbes};
            
            % Element numbers of those boundary elements that have only
            % Dirichlet boundary nodes. (Not needed?)
            this.Dir_boundary_elems = [];
            
            % Node numbers of all non-dirichlet-boundary nodes
            this.All_nondir_nodes = setdiff(1:this.N_node,All_Dir_boundary_nodes);
            
            % Map: node -> coordinate
            this.Node2coord = zeros(2,this.N_node);
            this.Node2coord(1,:) = reshape(this.Geo_node{1}',[1,numel(this.Geo_node{1})]);
            this.Node2coord(2,:) = reshape(this.Geo_node{2}',[1,numel(this.Geo_node{2})]);
            
            % Element number -> parameter
            this.Elem2param = PrimaryModel.PDEParameters.evaluateParameterField(...
                reshape(this.Geo_elem{1}',[1,numel(this.Geo_elem{1})]),...
                reshape(this.Geo_elem{2}',[1,numel(this.Geo_elem{2})]));
            
            % Map: node -> equation
            this.Node2eq = 1:(this.Nx+1)*(this.Ny+1);
        end