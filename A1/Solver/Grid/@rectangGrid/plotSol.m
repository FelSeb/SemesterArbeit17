function [] = plotSol(this,PostSol, Display_gradient)

figure;
Tgrid = PostSol.T(this.Topo_node');

surf(this.Geo_node{1}, this.Geo_node{2}, Tgrid);
xlabel('x')
ylabel('y')
zlabel('T')


if(Display_gradient && isfield(PostSol,'dT'))
    dT = PostSol.dT;
    loc = PostSol.loc;
    
    % Reshape to grid
    dTdxgrid = reshape(dT(1,this.Topo_elem),size(this.Topo_elem))';
    dTdygrid = reshape(dT(2,this.Topo_elem),size(this.Topo_elem))';
    
    
    grad_grid_x = reshape(loc(1,this.Topo_elem),size(this.Topo_elem))';
    grad_grid_y = reshape(loc(2,this.Topo_elem),size(this.Topo_elem))';
    
    figure;
    surf(grad_grid_x, grad_grid_y, dTdxgrid);
    xlabel('x')
    ylabel('y')
    zlabel('dTdx')
    
    figure;
    surf(grad_grid_x, grad_grid_y, dTdygrid);
    xlabel('x')
    ylabel('y')
    zlabel('dTdy')
end
end