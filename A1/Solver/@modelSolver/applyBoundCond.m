function [ this] = applyBoundCond( this, varargin )
% This function applies Dirichlet and v. Neumann boundary conditions to the
% edges of a RECTANGULAR domain. The right hand side of the linear system
% of equations is composed.

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

nel = this.Grid.N_elem;

% Compose Rhs
this.Rhs = zeros(length(this.C),1);
elem_vecs = zeros(this.Grid.FE.Nodes_per_elem,nel);

% In case the derivatives of the finite element force vector with respect
% to the material parameters are to be computed aswell...
if(compute_grad)
    loc_Rhs_grad = zeros(length(this.C), this.Grid.N_elem * nparam);
    elem_vecs_grad = zeros(this.Grid.FE.Nodes_per_elem,nparam,nel);
end


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
        
        % determine indices of matrix entries
        
        loc_ind = this.Grid.Neu_boundary_elems{2}{edge}(:,bound_elem);
        
        % integral (2D Gauss quadrature)
        for p = 1:length(this.Grid.FE.NumInt1D.in)
            
            % Compute Jacobian matrix and its determinant
            % Assuming (bi-)linear basis functions here!!!
            J = norm(xu-xl)/2;
            detJ = J;
            
            N_A  = this.N_all1D(p, :)';
            
            wu = (this.Grid.FE.NumInt1D.in(p)+1)/2;
            wl = 1-wu;
            x = xl * wl + xu * wu;
            
            dTdx = this.PrimaryModel.Boundary_conds.Bound_cond_funs{2,edge}(x(1),x(2));
            
            elem_vecs(loc_ind,elem) = elem_vecs(loc_ind,elem) + ...
                this.Grid.FE.NumInt1D.weights(p) * N_A * (lambda * dTdx)'...
                * this.PrimaryModel.Domain.All_edges_normvec(:,edge) * detJ;
            
            
            % Compute derivative
            if(compute_grad)
                inds = (elem-1)*nparam+1 : (elem-1)*nparam+nparam;
                Xe = this.Grid.X(inds);
                elem_vecs_grad(loc_ind,:,elem) = elem_vecs_grad(loc_ind,:,elem) + ...
                    this.Grid.FE.NumInt1D.weights(p) * detJ * N_A* ...
                    getElemForceGrad(Xe, dTdx, this.PrimaryModel.Domain.All_edges_normvec(:,edge), this.Param_trafo);
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

% Apply Dirichlet BCs

% Method that returns dirichlet boundary values and other things in
% a convenient format
[dirs, dir_vals, nondir] = this.formatDirBoundInfo;

% Modifying the system of equations
this.C = this.C(nondir,:);
this.Rhs = this.Rhs(nondir) - this.C(:,dirs) * dir_vals;
this.C = this.C(:,nondir);



if(compute_grad)
    for elem  = this.Grid.All_boundary_elems
        lin = this.Grid.Node2eq(this.Grid.Elem2node(:, elem));
        cols = (elem-1)*nparam+1:(elem-1)*nparam+nparam;
        loc_Rhs_grad(lin,cols) = loc_Rhs_grad(lin,cols) + elem_vecs_grad(:,:,elem);
    end
    
    % convert to cell array (this must be done since each column of the
    % Rhs_grad will be shortened by as many entries as there are dirichlet
    % boundary nodes):
    loc_Rhs_grad = mat2cell(loc_Rhs_grad,size(loc_Rhs_grad,1),ones(1,this.Grid.N_elem * nparam));
    
    % Apply Dirichlet BCs
    for nx = 1:numel(this.C_grad)
        % Modifying the system of equations
        this.C_grad{nx} = this.C_grad{nx}(nondir,:);
        loc_Rhs_grad{nx} = loc_Rhs_grad{nx}(nondir) - this.C_grad{nx}(:,dirs) * dir_vals;
        this.C_grad{nx} = this.C_grad{nx}(:,nondir);
    end
    
    % Converting back Rhs_grad from cell to conventional array:
    this.Rhs_grad = cell2mat(loc_Rhs_grad);
end



end






function dfdX = getElemForceGrad(Xe, dTdx, normal, paramTrafo)

nparam = length(Xe);
if(nparam == 1)
    dfdl(1) = dTdx(1)*normal(1) + dTdx(2)*normal(2);
elseif(nparam == 2)
    dfdl(1) = dTdx(1)*normal(1);
    dfdl(2) = dTdx(2)*normal(2);
else
    dfdl(1) = dTdx(1)*normal(1);
    dfdl(2) = dTdx(2)*normal(2);
    dfdl(3) = dTdx(2)*normal(1) + dTdx(1)*normal(2);
end


dldX = paramTrafo.getDldXe(Xe);

dfdX = zeros(1,nparam);
for ix = 1:nparam
    for il = 1:nparam
        dfdX(ix) = dfdX(ix) + dfdl(il) * dldX(il,ix);
    end
end

end


