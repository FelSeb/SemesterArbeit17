classdef DataPoint
    %DATAPOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % input (in Elem2param format)
        x;
        
        % output (in PostSol.T format)
        y;
        
        %
        Info;
    end
    
    methods
        % Constructor
        function this = DataPoint(solver)
            this.Info = struct();
            this.Info.GridInfo = solver.Grid.Info;
            this.Info.DomainInfo = solver.PrimaryModel.Domain.Info;
            this.Info.ParamFieldInfo = solver.PrimaryModel.PDEParameters.Info;
            this.Info.BoundaryInfo = solver.PrimaryModel.Boundary_conds.Info;
        
            this.x = solver.Grid.Elem2param;
            this.y = solver.PostSol.T;
        end   
    end
    
end

