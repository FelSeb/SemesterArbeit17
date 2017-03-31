classdef PropDensity
    %PROPOSER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        Prop_step_length;
        
    end
    
    methods (Abstract)
        pval = proppdf(this, X_new, X_old, varargin);
        
        logpval = logproppdf(this, X_new, X_old, varargin);
        
        smpl = proprnd(this, X_old, varargin);
    end
    
end

