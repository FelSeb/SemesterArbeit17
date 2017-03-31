function [dirs, dir_vals, nondir] = formatDirBoundInfo(this)

% All equations corresponding to nondir nodes
nondir = this.Grid.Node2eq(this.Grid.All_nondir_nodes);

% Number of dirichlet nodes per edge
N_dirs = zeros(1,this.PrimaryModel.Domain.N_edges);
for edge= 1:this.PrimaryModel.Domain.N_edges
    N_dirs(edge) = length(this.Grid.Dir_boundary_nodes{edge});
end

% Composing the voector of prescribed boundary values
dir_vals = zeros(sum(N_dirs),1);
dirs = zeros(sum(N_dirs),1);
from = 1;
for edge= 1:this.PrimaryModel.Domain.N_edges
    % Boundary edge
    to = from + N_dirs(edge) - 1;
    
    dir_coords = this.Grid.Node2coord(:,this.Grid.Dir_boundary_nodes{edge});
    dir_vals(from:to) = this.PrimaryModel.Boundary_conds.Bound_cond_funs{1,edge}(dir_coords(1),dir_coords(2));
    dirs(from:to) = this.Grid.Dir_boundary_nodes{edge};
    
    from = to + 1;
    if(from > sum(N_dirs))
       break; 
    end
end



end

