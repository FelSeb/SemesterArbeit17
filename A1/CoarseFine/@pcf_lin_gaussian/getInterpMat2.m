function [ this ] = getInterpMat2( this,isoshape )

% Bilinear interpolation (arbitrary quadrilateral !!!):

if(isoshape) % if all CG elements have the same shape
    % Coordinates of the quadrilateral
    x = [this.Coarse_grid{1}(1,1), this.Coarse_grid{2}(1,1);
        this.Coarse_grid{1}(1,2), this.Coarse_grid{2}(1,2);
        this.Coarse_grid{1}(2,2), this.Coarse_grid{2}(2,2);
        this.Coarse_grid{1}(2,1), this.Coarse_grid{2}(2,1); ]';
    
    
    % Coordinate for interpolation
    x_star(1,:) = reshape(this.Fine_grid{1},1,numel(this.Fine_grid{1}));
    x_star(2,:) = reshape(this.Fine_grid{2},1,numel(this.Fine_grid{2}));
    
    
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
    
    N = zeros(4,length(x_star(1,:)));
    for i = 1:length(x_star(1,:))
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
        
        
        N1 = @(xi,eta) 0.25 * (1-xi-eta+xi*eta);
        N2 = @(xi,eta) 0.25 * (1-xi+eta-xi*eta);
        N3 = @(xi,eta) 0.25 * (1+xi+eta+xi*eta);
        N4 = @(xi,eta) 0.25 * (1+xi-eta-xi*eta);
        
        if(xi(1) > 1 ||xi(1) < -1 || xi(2) > 1 || xi(2) < -1)
            N(:,i) = [0;0;0;0];
        else
            N(:,i) = [N1(xi(1),xi(2)); N2(xi(1),xi(2)); N3(xi(1),xi(2)); N4(xi(1),xi(2))];
        end
        
    end
    
    this.Interpolation_matrix = N';
else
    
    
end
end

