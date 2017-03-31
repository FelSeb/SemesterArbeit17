function [ this] = applyBoundCond2( this, varargin )
warning('verwende eher applyBoundCond3')


% This function applies Dirichlet and v. Neumann boundary conditions to the
% edges of a RECTANGULAR domain.

% Apply v. Neumann BCs

E2Pdim = size(this.Grid.Elem2param);
if(length(E2Pdim) == 3)
    nparam = 3;
elseif(E2Pdim(1) == 2)
    nparam = 2;
else
    nparam = 1;
end

if(numel(varargin)>0)
    compute_grad = varargin{1};
else
    compute_grad = 0;
end

% Compose Rhs
this.Rhs = zeros(length(this.C),1);
if(compute_grad)
    this.Rhs_grad = sparse(length(this.C), this.Grid.N_elem * nparam);
end

% To Do: vercorize operations


nel = this.Grid.N_elem;
elem_vecs = zeros(4,nel);
elem_vecs_grad = zeros(4,nparam,nel);

for edge= 1:this.PrimaryModel.Domain.N_edges
    % Boundary edge
    
    bound_elem = 1;
    for elem  = this.Grid.Neu_boundary_elems{1}{edge}
        
        % the coordinates of the element nodes
        xe = this.Grid.Node2coord(:,this.Grid.Elem2node(:,elem));
        
        xl = xe(:,this.Grid.Neu_boundary_elems{2}{edge}(1,bound_elem));
        xu = xe(:,this.Grid.Neu_boundary_elems{2}{edge}(2,bound_elem));
        
        
        % get local material parameter (distiguish 3 different formats)
        if(nparam == 3)
            lambda = this.Grid.Elem2param(:,:,elem);
        elseif(nparam == 2)
            lambda = diag(this.Grid.Elem2param(:,elem));
        else
            lambda = this.Grid.Elem2param(elem);
        end
        
        % Compose element matrix
        for nodA = [1,2] % Test
            
            % determine indices of matrix entries
            lin = this.Grid.Node2eq(this.Grid.Elem2node( ...
                this.Grid.Neu_boundary_elems{2}{edge}(nodA,bound_elem), elem ));
            
            loc_ind = this.Grid.Neu_boundary_elems{2}{edge}(nodA,bound_elem);
            
            % integral (2D Gauss quadrature)
            for p = 1:length(this.Grid.FE.NumInt1D.in)
                
                % Compute Jacobian matrix and its determinant
                J = norm(xu-xl)/2;
                detJ = J;
                
                N_A  = this.N_all1D(p, nodA);
                
                wu = (this.Grid.FE.NumInt1D.in(p)+1)/2;
                wl = 1-wu;
                x = xl * wl + xu * wu;
                
                dTdx = this.PrimaryModel.Boundary_conds.Bound_cond_funs{2,edge}(x(1),x(2));
                
                elem_vecs(loc_ind,elem) = elem_vecs(loc_ind,elem) + ...
                    this.Grid.FE.NumInt1D.weights(p) * N_A * (lambda * dTdx)'...
                    * this.PrimaryModel.Domain.All_edges_normvec(:,edge) * detJ;
                
                
                % Compute derivative
                if(compute_grad)
                    %error('this is wrong: remeber that the derivative is with respect to all individual elements in X')
                    elem_vecs_grad(loc_ind,1,elem) = elem_vecs_grad(loc_ind,1,elem) + ...
                        this.Grid.FE.NumInt1D.weights(p) * N_A * lambda(1,1) * dTdx(1)...
                        * this.PrimaryModel.Domain.All_edges_normvec(1,edge) * detJ;
                    if(nparam > 1)
                        elem_vecs_grad(loc_ind,2,elem)  = elem_vecs_grad(loc_ind,2,elem) + ...
                            this.Grid.FE.NumInt1D.weights(p) * N_A * lambda(2,2) * dTdx(2)...
                            * this.PrimaryModel.Domain.All_edges_normvec(2,edge) * detJ;
                    end
                    if(nparam == 3)
                        elem_vecs_grad(loc_ind,3,elem) = elem_vecs_grad(loc_ind,3,elem) + ...
                            this.Grid.FE.NumInt1D.weights(p) * N_A * lambda(1,2) * ...
                            (dTdx(2) * this.PrimaryModel.Domain.All_edges_normvec(1,edge) + ...
                            dTdx(1) * this.PrimaryModel.Domain.All_edges_normvec(2,edge))...
                            * detJ;
                    end
                end
                
                
            end
        end
        bound_elem = bound_elem+1;
    end
end

% Compose global vector
for elem  = this.Grid.All_boundary_elems
    lin = this.Grid.Node2eq(this.Grid.Elem2node(:, elem));
    this.Rhs(lin) = this.Rhs(lin) + elem_vecs(:, elem);
end

if(compute_grad)
    for elem  = this.Grid.All_boundary_elems
        lin = this.Grid.Node2eq(this.Grid.Elem2node(:, elem));
        cols = (elem-1)*nparam+1:(elem-1)*nparam+nparam;
        this.Rhs_grad(lin,cols) = this.Rhs_grad(lin,cols) + elem_vecs_grad(:,:,elem);
    end
    % convert to cell array:
    this.Rhs_grad = mat2cell(this.Rhs_grad,length(this.C),ones(1,this.Grid.N_elem * nparam));
end


% Apply Dirichlet BCs

% Method that returns dirichlet boundary values and other things in
% a convenient format
[dirs, dir_vals, nondir] = this.formatDirBoundInfo;

% Modifying the system of equations
this.C = this.C(nondir,:);
this.Rhs = this.Rhs(nondir) - this.C(:,dirs) * dir_vals;
this.C = this.C(:,nondir);

if(compute_grad)
    for nx = 1:numel(this.C_grad)
        this.C_grad{nx} = this.C_grad{nx}(nondir,:);
        % Modifying the system of equations
        this.Rhs_grad{nx} = this.Rhs_grad{nx}(nondir) - this.C_grad{nx}(:,dirs) * dir_vals;
        this.C_grad{nx} = this.C_grad{nx}(:,nondir);
    end
end

end

