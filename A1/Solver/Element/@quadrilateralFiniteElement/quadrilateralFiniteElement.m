classdef quadrilateralFiniteElement < finiteElement
    %QUADRIALTERALFINITEELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nodes_per_elem;
        Nodes_per_elem1D;
        
        % Shape- /test- functions (and derivatives)
        SF;
        SF1D;
        
        % Numerical integration
        NumInt;
        NumInt1D;
    end
    
    methods
        % Constructor
        function this = quadrilateralFiniteElement(sf_string, numint_string)
            if(strcmp(sf_string,'sf_lin'))
                this.SF = @( xi, eta , nf) this.SF_lin( xi, eta , nf);
                this.SF1D = @( xi, nf) this.SF_lin1D( xi, nf);
                this.Nodes_per_elem = 4;
                this.Nodes_per_elem1D = 2;
            else
                error('The selected shape-/ test- function does not exist')
            end
            
            if(strcmp(numint_string,'GaussQuad3'))
                [this.NumInt, this.NumInt1D] = this.setUpNumIntGauss3;
            elseif(strcmp(numint_string,'GaussQuad2'))
                [this.NumInt, this.NumInt1D] = this.setUpNumIntGauss2;
            else
                error('The selected num. int.- scheme is not available')
            end
        end
        
    end
    
    methods (Static)
        
        % Linear shape and test functions
        [ val, dval ] = SF_lin( xi, eta , nf);
        [val, dval] =  SF_lin1D( xi, nf);
        
        % Set up 1- and 2D numerical integration schemes
        function [NumInt, NumInt1D] = setUpNumIntGauss3()
            %Gauss quadrature input points
            [NumInt.in{1}, NumInt.in{2}] = meshgrid([-sqrt(0.6),0,sqrt(0.6)]);
            % Gauss quadrature weights
            NumInt.weights = [25,40,25; 40,64,40; 25,40,25]/81;
            
            % 1D Numerical intergation
            NumInt1D.in = [-sqrt(3/5),0,sqrt(3/5)]; % Gauss quadrature points
            NumInt1D.weights = [5/9,8/9,5/9];% Gauss quadrature weights
        end
        
                % Set up 1- and 2D numerical integration schemes
        function [NumInt, NumInt1D] = setUpNumIntGauss2()
            %Gauss quadrature input points
            [NumInt.in{1}, NumInt.in{2}] = meshgrid([-sqrt(1/3),sqrt(1/3)]);
            % Gauss quadrature weights
            NumInt.weights = [1, 1; 1, 1];
            
            % 1D Numerical intergation
            NumInt1D.in = [-sqrt(1/3),sqrt(1/3)]; % Gauss quadrature points
            NumInt1D.weights = [1,1];% Gauss quadrature weights
        end
        
        
        function [ e2n] = GetElem2node(grid)
            % Get the element2node map for quadrilateral elements
            % Usage:
            % [ e2n ] = GetElem2nod( Grid )
            
            if(isa(grid,'rectangGrid'))
                % 1 o---o 2
                %   |   |
                % 4 o---o 3
                
                e2n = zeros(numel(grid.Topo_elem),4);
                for i = 1:grid.Nx
                    for j = 1:grid.Ny
                        e2n(grid.Topo_elem(i,j),1) = grid.Topo_node(i,j);
                        e2n(grid.Topo_elem(i,j),2) = grid.Topo_node(i+1,j);
                        e2n(grid.Topo_elem(i,j),3) = grid.Topo_node(i+1,j+1);
                        e2n(grid.Topo_elem(i,j),4) = grid.Topo_node(i,j+1);
                    end
                end
                e2n = e2n';
            else
                error('Unknown grid type')
            end
        end
        
       
        
    end
    
end

