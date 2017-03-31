classdef pdeParameters
    %PDE_PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Array containing the thermal conductivities of all materials of
        % the compound. Cell array in case of anisotropy.
        Lambdas;
        
        % Amount of different materials in the compound
        Material_multiplicity;
        
        % String specifying according to which rule the parameter field is
        % to be generated
        Parameter_field;
    end
    
    methods
        % Constructor
        function this = pdeParameters(lambdas, cps, rhos, param_field_string)
            this.Lambdas = lambdas;
            
            if(length(lambdas) == length(cps) && length(rhos) == length(cps))
                this.Material_multiplicity = length(lambdas);
            else
                error('labdas, cps and rhos must be of the same length')
            end
            
            this.Parameter_field = param_field_string;
        end
        
        % Set parameter field type
        function this = setParamField(this, param_field_string)
            this.Parameter_field = param_field_string;
        end
        
        % Evalate the parameter field
        function Lambda_grid = evaluateParameterField(this, x_grid,y_grid)
            if(strcmp(this.Parameter_field,'deterministic1'))
                Lambda_grid = deterministicField1(this, x_grid,y_grid);
            elseif(strcmp(this.Parameter_field,'deterministic2'))
                Lambda_grid = deterministicField2(this, x_grid,y_grid);
            elseif(strcmp(this.Parameter_field,'random1'))
                Lambda_grid = randomField1(this, x_grid,y_grid);
            end
        end
        
        
        % Declaring various parameter field functions
        Lambda_grid = deterministicField1(this, x_grid,y_grid);
        
        Lambda_grid = deterministicField2(this, x_grid,y_grid);
        
        Lambda_grid = randomField1(this, x_grid,y_grid);
        
        % Plotting the parameter field for a given coordinate input in
        % meshgrid format.
        function plotField(this,xg,yg)
           params = this.evaluateParameterField(...
               reshape(xg, [1,numel(xg)]),...
               reshape(yg, [1,numel(yg)]));
           colorp(xg,yg,reshape(params,size(xg)));
        end
        
    end
    
    
    
end

