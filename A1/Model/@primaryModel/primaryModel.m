classdef primaryModel
    %PRIAMARYMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Domain; % Defines the domains geometry and boundary
        
        PDEParameters; % Defines the (random) spatial parameter field
        
        Boundary_conds; % defines the bondary conditions
    end
    
    methods
        function this = primaryModel(domain, pde_parameters, bounary_conds)
            this.Domain = domain;
            this.PDEParameters = pde_parameters;
            this.Boundary_conds = bounary_conds; 
        end
    end
    
end

