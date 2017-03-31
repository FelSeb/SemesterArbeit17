function [ this ] = solveSystem( this, varargin )

this.Sol = this.C\this.Rhs;

if(numel(varargin)>0)
    compute_sol_grad = varargin{1};
else
    compute_sol_grad = 0;
end

% Not needed:
if(compute_sol_grad)
    % determine Elem2param format
    E2Pdim = size(this.Grid.Elem2param);
    if(length(E2Pdim) == 3)
        nparam = 3;
    elseif(E2Pdim(1) == 2)
        nparam = 2;
    else
        nparam = 1;
    end
    
    % Bad Idea: 
    for np = 1:nparam * this.Grid.N_elem
        this.Sol_grad(:,np) = this.C \ (-this.C_grad{np} * this.Sol + this.Rhs_grad(:,np));
    end
    
end


end

