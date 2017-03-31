classdef (Abstract) finiteElement 
    %FINITEELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        Nodes_per_elem;
        
        % Shape- /test- functions (and derivatives)
        SF;
        
        % Numerical integration
        NumInt;
        NumInt1D;
        
    end
    
    
end

