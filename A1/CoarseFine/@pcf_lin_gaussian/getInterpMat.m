function [ this ] = getInterpMat( this,isoshape )
% Using a bilinear function for the interpolation of a quadrialteral can be
% achieved by solving a quadratic equation for each of the interpolation
% points. The general solution can be obtained by hand and then be used to
% compute the entries of the interpolation matrix.

% f

% Bilinear interpolation (arbitrary quadrilateral !!!):
if(~isa(this.Coarse_grid.FE,'quadrilateralFiniteElement'))
    error('getInterpMat can only be applied to quadrilateral coarse elements')
end

N1 = @(xi,eta) 0.25 * (1-xi-eta+xi*eta);
N2 = @(xi,eta) 0.25 * (1+xi-eta-xi*eta);
N3 = @(xi,eta) 0.25 * (1+xi+eta+xi*eta);
N4 = @(xi,eta) 0.25 * (1-xi+eta-xi*eta);

N = zeros(this.Fine_grid.N_node, this.Coarse_grid.N_node);

% In order to avoid that some rows of the interpolation matrix are
% computed multiple times rows numbers that are already computed are
% saved in the array covered_nodes.
covered_nodes = [];
for ce = 1:this.Coarse_grid.N_elem
    
    % Coordinates of the coarse element
    coarse_nodes = this.Coarse_grid.Elem2node(:,ce);
    x = this.Coarse_grid.Node2coord(:,coarse_nodes);
    
    % Some precomputation
    A = 0.25 * (x(1,1)+x(1,2)+x(1,3)+x(1,4));
    B = 0.25 * (-x(1,1)-x(1,2)+x(1,3)+x(1,4));
    C = 0.25 * (-x(1,1)+x(1,2)+x(1,3)-x(1,4));
    D = 0.25 * (x(1,1)-x(1,2)+x(1,3)-x(1,4));
    E = 0.25 * (x(2,1)+x(2,2)+x(2,3)+x(2,4));
    F = 0.25 * (-x(2,1)-x(2,2)+x(2,3)+x(2,4));
    G = 0.25 * (-x(2,1)+x(2,2)+x(2,3)-x(2,4));
    H = 0.25 * (x(2,1)-x(2,2)+x(2,3)-x(2,4));

    
    % coefficient a of quadratic formula
    a = (F*D -H*B);
    
    % Fine Grained nodes
    all_fine_elems = this.C2F{ce};
    all_fine_nodes = unique(reshape( ...
        this.Fine_grid.Elem2node(:,all_fine_elems),...
        [1,length(all_fine_elems)*4]...
        ));
    fine_nodes = setdiff(all_fine_nodes, covered_nodes);
    covered_nodes = [covered_nodes,fine_nodes];
    
    % Coordinates of the fine grained nodes
    x_star = this.Fine_grid.Node2coord(:,fine_nodes);
    
    
    for i = 1:length(fine_nodes)
        % coefficient b of quadratic formula
        b = (-x_star(2,i)*D +E*D +F*C -G*B -H*A + H*x_star(1,i));
        % coefficient c of quadratic formula
        c = (E -x_star(2,i))*C - (A - x_star(1,i))*G;
        
        if(a~=0)
            xi(2) = (-b + sqrt(b^2-4*a*c))/(2*a); % OR: (-b - sqrt(b^2-4*a*c))/(2*a);
            xi(1) = (x_star(1,i)-A-B*xi(2))/(C+D*xi(2));
        else
            xi(2) = -c/b;
            xi(1) = (x_star(1,i)-A-B*xi(2))/(C+D*xi(2));
        end
        
        if(xi(1) > 1+1e-12 ||xi(1) < -1-1e-12 || xi(2) > 1+1e-12 || xi(2) < -1-1e-12)
            N(fine_nodes(i),coarse_nodes) = [0,0,0,0];
        else
            N(fine_nodes(i),coarse_nodes) = [N1(xi(1),xi(2)), N2(xi(1),xi(2)), N3(xi(1),xi(2)), N4(xi(1),xi(2))];
        end  
    end
end

this.Interpolation_matrix = N;
end

