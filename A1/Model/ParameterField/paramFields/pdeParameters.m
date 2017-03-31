classdef pdeParameters
    %PDE_PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Array containing the thermal conductivities of all materials of
        % the compound. Cell array in case of anisotropy.
        Lambdas;
        
        % Amount of different materials in the compound
        Material_multiplicity;
    end
    
    methods
        % Constructor
        function this = pdeParameters(lambdas, cps, rhos)
            this.Lambdas = lambdas;
            
            if(length(lambdas) == length(cps) && length(rhos) == length(cps))
                this.Material_multiplicity = length(lambdas);
            else
                error('labdas, cps and rhos must be of the same length')
            end
            
        end
        
        % Plotting the parameter field for a given coordinate input in
        % meshgrid format.
        function params = plotField(this,xg,yg)
           params = this.evaluateParameterField(...
               reshape(xg, [1,numel(xg)]),...
               reshape(yg, [1,numel(yg)]));
           pcolor(xg,yg,reshape(params,size(xg)));
        end
        
    end
    
    methods (Abstract)
        % Evalate the parameter field
        Lambdas = evaluateParameterField(this, x_grid,y_grid)
    end
end