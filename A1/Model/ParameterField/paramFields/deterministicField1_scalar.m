classdef deterministicField1_scalar < pdeParameters
    %DETERMINITICFIELD1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Info;
    end
    
    methods
        % Constructor
        function this = deterministicField1_scalar(lambdas, cps, rhos)
            this = this@pdeParameters(lambdas, cps, rhos);
            
            this.Info = struct;
            this.Info.Type = class(this);
        end
        
        function [ lambdas ] = evaluateParameterField(this, x_coords,y_coords )
            % Uniform, isotropic parameter.
            
            [row,col] = size(x_coords);
            
            if(row > 1)
                error('coordinate inputs to deterministicField1 have wrong format')
            end
            
            lambdas = ones(length(x_coords),1) * this.Lambdas{1};
            
        end
        
    end
    
end

