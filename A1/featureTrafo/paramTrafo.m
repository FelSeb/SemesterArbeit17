classdef paramTrafo
    %PARAMTRAFO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Abstract)
        
        dldX = getDldXe(Xe)

        
        X = getXfromElem2param(Elem2param)
  
        
        lambda = getlformXe(Xe)

        
        
        e2p = getElem2ParamfromX(X,nelem)

        
        
    end
    
end

