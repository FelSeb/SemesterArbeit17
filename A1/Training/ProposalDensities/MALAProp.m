classdef MALAProp < PropDensity
    %MALAPROP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Proposal variance
        Prop_step_length;
        
        Step_length = 1;
    end
    
    methods
        
        % Constructor (empty)
        function this = MALAProp(prop_var)
            this.Prop_step_length = prop_var;
        end
        
        
        function pval = proppdf(this, X_new, X_old, varargin)
            pval = 0;
            error('not implemented');
        end
        
        
        function logp = logproppdf(this, x, y, varargin)
            if(numel(varargin) > 0)
                info =  varargin{1};
            else
                info = cell(1,2);
                info(:) = {0};
                warning(['MALALogProp not provided with gradient of logpdf.',...
                    ' This is not problematic if MALALogProp is only tested in the beginning of MYmhsample'])
            end
            
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
            
            
            sigma2 = this.Prop_step_length * this.Step_length;
            mean = y' + sigma2/2 * info{1,2};
            % cova = eye(length(y))*sigma2;
            % logp = log(mvnpdf(x',mean,cova));
            
            k = length(mean);
            vec = (x-mean');
            logp = -k/2*log(2*pi) -0.5*k*log(sigma2) -0.5*(vec' *vec)/sigma2;
        end
        
        function smpl = proprnd(this, y, varargin)
            % make sure y is a column vector
            if(length(y) == numel(y))
                y = reshape(y,numel(y),1);
            else
                error('y has wrong format')
            end
            
            if(numel(varargin) > 0)
                info =  varargin{1};
            else
                info = cell(1,2);
                info(:) = {0};
                warning('MALALogProp not provided with gradient of logpdf.')
            end

            sigma2 = this.Prop_step_length * this.Step_length;
            mean = y' + sigma2/2 * info{1,2};
            cova = eye(length(y))*sigma2;
            
            smpl = mvnrnd(mean,cova);
        end
        
        
    end
    
end

