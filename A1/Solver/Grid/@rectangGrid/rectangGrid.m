classdef rectangGrid < simplegrid
    % This is the class of rectangular grids. It is a simple grid structure
    % that allows certain simplifications in the algorithm. It can only be
    % used for rectangular Domains.
    
    properties
        Topo_elem; % Topological grid of elements
        Topo_node; % Topological grid of nodes
        Geo_elem; % Geometrical positions of all element centers
        Geo_node; % Geometrical positions of all nodes
        
        
        % Number of elements in X direction
        Nx;
        % Number of elements in Y direction
        Ny;
        
        % Information concerning the grid
        Info;
    end
    
    methods
        
        % Constructor
        function this = rectangGrid(nx,ny,PrimaryModel,fe_type_string,sf_string,numint_string)
            if(nx > 0 && ny > 0)
                this.Topo_elem = reshape(1:nx*ny, nx, ny);
                this.Topo_node = reshape(1:(nx+1)*(ny+1),[nx+1, ny+1]);
                
                xmin = PrimaryModel.Domain.X_min;
                xmax = PrimaryModel.Domain.X_max;
                ymin = PrimaryModel.Domain.Y_min;
                ymax = PrimaryModel.Domain.Y_max;
                dx = (xmax-xmin)/nx;
                dy = (ymax-ymin)/ny;
                
                this.Geo_node{1} = zeros(nx+1, ny+1);
                this.Geo_node{2} = zeros(nx+1, ny+1);
                [this.Geo_node{1},this.Geo_node{2}] = ...
                    meshgrid(   linspace(xmin,xmax,nx+1), linspace(ymin,ymax,ny+1));
                
                this.Geo_elem{1} = zeros(nx, ny);
                this.Geo_elem{2} = zeros(nx, ny);
                [this.Geo_elem{1},this.Geo_elem{2}] = ...
                    meshgrid(   linspace(xmin+dx/2,xmax-dx/2,nx), linspace(ymin+dy/2,ymax-dy/2,ny));
                
                this.N_elem = nx * ny;
                this.N_node = (nx+1) * (ny+1);
                this.Nx = nx;
                this.Ny = ny;
                
                if(strcmp(fe_type_string,'quadrilateral'))
                    this.FE = quadrilateralFiniteElement(sf_string, numint_string);
                else
                    error('The element type does not exist')
                end
                
                % Maps:
                this = this.makeMaps(PrimaryModel);
                
                % Store some information
                this.Info = struct;
                this.Info.Nx = this.Nx;
                this.Info.Ny = this.Ny;
                this.Info.FE_type = class(this.FE);
                this.Info.ShapeFun_type = sf_string;
                this.Info.NumIntScheme = numint_string;
            else
                
            end
        end
        
        % The makeMaps functions constructs the following maps:
        % Elem2node
        % Node2coord
        % Elem2param
        % Node2eq
        % All_boundary_nodes
        % All_boundary_elems
        % Dir_boundary_nodes
        % Neu_boundary_nodes
        % Neu_boundary_elems
        % Dir_boundary_elems
        % All_nondir_nodes
        this = makeMaps(this, PrimaryModel);
        
        
        % If the material parameter field is random this function can be
        % used to obtain a new random sample.
        function this = remakeElem2Param(this, PrimaryModel)
            % Element number -> parameter
            this.Elem2param = PrimaryModel.PDEParameters.evaluateParameterField(...
                reshape(this.Geo_elem{1}',[1,numel(this.Geo_elem{1})]),...
                reshape(this.Geo_elem{2}',[1,numel(this.Geo_elem{2})]));
        end
        
        
        % The makeC2F function constructs a map that which fine elements
        % belong to which coarse element
        c2f = makeC2F(this, grid2)
        
        % The plotSol function creates a plot for a given postprocessed
        % solution
        [] = plotSol(this, PostSol, Display_gradient)
        
        % Construct certain maps corresponding to a coarse grained element
        % which comprises the elements contained in local_elem_numbers.
        function [sub_maps,ce_node_coords] = getSubMaps(this, local_elem_numbers)
            isSub = ismember(this.Topo_elem,local_elem_numbers);
            [a,b] = find(isSub == 1);
            au = unique(a);
            bu = unique(b);
            nx = length(au);
            ny = length(bu);
            
            sub_maps.Topo_elem = reshape(1:nx*ny, nx, ny);
            sub_maps.Geo_elem{1} = this.Geo_elem{1}(au,bu);
            sub_maps.Geo_elem{2} = this.Geo_elem{2}(au,bu);
            sub_maps.Nx = nx;
            sub_maps.Ny = ny;
            
            if(nargout > 1)
                % elements of the fine grid that have nodes in common with the
                % coarse element.
                cn_elems(1) = this.Topo_elem(au(1),bu(1));
                cn_elems(2) = this.Topo_elem(au(end),bu(1));
                cn_elems(3) = this.Topo_elem(au(end),bu(end));
                cn_elems(4) = this.Topo_elem(au(1),bu(end));
                
                ce_node_coords(:,1) = this.Node2coord(:,this.Elem2node(1,cn_elems(1)));
                ce_node_coords(:,2) = this.Node2coord(:,this.Elem2node(2,cn_elems(2)));
                ce_node_coords(:,3) = this.Node2coord(:,this.Elem2node(3,cn_elems(3)));
                ce_node_coords(:,4) = this.Node2coord(:,this.Elem2node(4,cn_elems(4)));
            end
            
            
        end
        
    end
end

