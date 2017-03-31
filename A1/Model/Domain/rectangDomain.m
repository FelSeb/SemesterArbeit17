classdef rectangDomain < polyangularDomain
    %RECTANGULARDOMAIN 
    
    properties
        
    % Information about the domain
    Info;
    end
    
    methods
        % Make a rectangular domain
        function this = rectangDomain(xmin,xmax,ymin,ymax)
            this.Edges = zeros(4,4);
            this.Edges(:,1) = [xmin,ymin,xmax,ymin]'; % a
            this.Edges(:,2) = [xmax,ymin,xmax,ymax]'; % b
            this.Edges(:,3) = [xmax,ymax,xmin,ymax]'; % c
            this.Edges(:,4) = [xmin,ymax,xmin,ymin]'; % d
            
            this.X_min = xmin;
            this.X_max = xmax;
            this.Y_min = ymin;
            this.Y_max = ymax;
            
            this.N_edges = 4;
            
            this = this.GetNormVecs;
            
            % set infos
            this.Info = struct;
            this.Info.X_min = xmin;
            this.Info.X_max = xmax;
            this.Info.Y_min = ymin;
            this.Info.Y_max = ymax;
        end
    end
    
    
end

