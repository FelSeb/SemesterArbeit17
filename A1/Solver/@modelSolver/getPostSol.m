function [this] = getPostSol( this, Compute_gradient, varargin )
% Usage:
% modelSolverObj = modelSolverObj.getPostSol( this, Compute_gradient, Format_to_grid )
% If the gradient of the solution at all element centers
% (i.e. for xi = eta = 0) is to be computed set Compute_gradient = 1.

%% Format temperature

% Method that returns dirichlet boundary values and other things in
% a convenient format
[dirs, dir_vals, inner] = this.formatDirBoundInfo;

Tt = zeros(numel(inner) + numel(dirs),1);

% Add Dirichlet BC node values
Tt(inner) = this.Sol;
Tt(dirs) = dir_vals;

this.PostSol.T = Tt;


%% Format temperature gradient

if(Compute_gradient)
    
    N = zeros(4,1);
    DN = zeros(4,2);
    
    for node = 1:4
        [N(node), DN(node,:)] = this.Grid.FE.SF(0,0, node);
    end
    
    dT = zeros(2,this.Grid.N_elem);
    loc = zeros(2,this.Grid.N_elem);
    
    parfor elem  = 1:this.Grid.N_elem
        % the coordinates of the element nodes
        xe = this.Grid.Node2coord(:,this.Grid.Elem2node(:,elem));
        loc(:,elem) = xe * N;
        
        J = xe * DN;
        detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1);
        invJ = 1/detJ * [J(2,2), -J(1,2); -J(2,1), J(1,1)];
        
        T_local = Tt(this.Grid.Node2eq(this.Grid.Elem2node(:,elem)));
        DN_middle = DN' * T_local;
        
        dT(:,elem) = invJ' * DN_middle;
    end
    
    this.PostSol.dT = dT;
    this.PostSol.loc = loc;
end


if(numel(varargin)>0)
    compute_Xgrad_post = varargin{1};
else
    compute_Xgrad_post = 0;
end


if(compute_Xgrad_post)
    % determine Elem2param format
    E2Pdim = size(this.Grid.Elem2param);
    if(length(E2Pdim) == 3)
        nparam = 3;
    elseif(E2Pdim(1) == 2)
        nparam = 2;
    else
        nparam = 1;
    end
    dTt = zeros(numel(inner) + numel(dirs),nparam * this.Grid.N_elem);
    % Add Dirichlet BC node values
    dTt(inner,:) = this.Sol_grad;
    dTt(dirs,:) = 0;
    
    this.PostSol.dYdX = dTt;
    
end


end