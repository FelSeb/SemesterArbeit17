classdef GaussianProp < PropDensity
    %GAUSSIANPROP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Proposal variance
        Prop_step_length;
    end
    
    methods
        
        % Constructor (empty)
        function this = GaussianProp(prop_var)
            this.Prop_step_length = prop_var;
        end
            
        
        function pval = proppdf(this, X_new, X_old, varargin)
            pval = 0;
            error('not implemented'); 
        end

        function [ logp ] = logproppdf(this, x,y, varargin)
            % Gaussian proposal distribution. Probability density of choosing x given
            % y. x and y are column vectors.
            
            % make sure x and y are column vectors
            if(length(x) == numel(x))
               x = reshape(x,numel(x),1);
            else
                error('x has wrong format')
            end
            if(length(y) == numel(y))
               y = reshape(y,numel(y),1);
            else
                error('y has wrong format')
            end
            
            k = length(x);
            vec = (x-y);
            logp = -k/2 * log(2*pi) - 0.5*k*log(this.Prop_step_length) - 0.5 * (vec'*vec) / this.Prop_step_length;
        end
        
        function [ smpl ] = proprnd(this, y, varargin )
            % make sure y is a column vector
            if(length(y) == numel(y))
                y = reshape(y,numel(y),1);
            else
                error('y has wrong format')
            end
            
            % Gaussian proposal distribution. Probability density of choosing x given
            % y. x and y are column vectors.
            smpl = mvnrnd(y',eye(length(y))*this.Prop_step_length);
        end
        
        
    end
    
end

