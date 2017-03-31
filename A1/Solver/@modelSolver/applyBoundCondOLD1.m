function [ this] = applyBoundCond( this )

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


% Compose Rhs
this.Rhs = zeros(length(this.C),1);

To Do: rewrite the following code and devide it into an element vector and an assembly part

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
                
                this.Rhs(lin) = this.Rhs(lin) + ...
                    this.Grid.FE.NumInt1D.weights(p) * N_A * (lambda * dTdx)'...
                    * this.PrimaryModel.Domain.All_edges_normvec(:,edge) * detJ;
                
            end
        end
        bound_elem = bound_elem+1;
    end
end

% Apply Dirichlet BCs

% Method that returns dirichlet boundary values and other things in
% a convenient format
[dirs, dir_vals, nondir] = this.formatDirBoundInfo;

% Modifying the system of equations
this.C = this.C(nondir,:);
this.Rhs = this.Rhs(nondir) - this.C(:,dirs) * dir_vals;
this.C = this.C(:,nondir);


end


