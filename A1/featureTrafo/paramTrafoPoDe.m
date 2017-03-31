classdef paramTrafoPoDe < paramTrafo
    % Paramtrization to ensure positive definiteness of the conductivity
    % parameter.
    % The corresponding parametrization of lambda through X is done as follows:
    % 1.)   lambda_1 = exp(2 * X_1)
    % 2.)   lambda_1 = exp(2 * X_1)
    %       lambda_2 = exp(2 * X_2)
    % 3.)   lambda_1 = exp(2 * X_1)
    %       lambda_2 = exp(2 * X_2) + X_3^2
    %       lambda_3 = exp(2 * X_1) * X_3
    
    properties
    end
    
    methods
        % Constructor
        function this = paramTrafoPoDe()
            
        end
        
    end
    
    
    
    methods(Static)
        function dldX = getDldXe(Xe)
            nparam = length(Xe);
            
            if(nparam == 1)
                dldX(1,1) = exp( 2* Xe(1)) * 2;
            elseif(nparam == 2);
                dldX(1,1) = exp( 2* Xe(1)) * 2;
                dldX(1,2) = 0;
                dldX(2,1) = 0;
                dldX(2,2) = exp( 2* Xe(2)) * 2;
            elseif(nparam == 3)
                dldX(1,1) = exp( 2* Xe(1)) * 2;
                dldX(1,2) = 0;
                dldX(1,3) = 0;
                dldX(2,1) = 0;
                dldX(2,2) = exp( 2* Xe(2)) * 2;
                dldX(2,3) = 2 * X(3);
                dldX(3,1) = exp(Xe(1))*Xe(3);
                dldX(3,2) = 0;
                dldX(3,3) = exp(Xe(1));
            end
            
        end
        
        function X = getXfromElem2param(Elem2param)
            E2Pdim = size(Elem2param);
            nelem = E2Pdim(end);
            if(length(E2Pdim) == 3)
                nparam = 3;
            elseif(E2Pdim(1) == 2)
                nparam = 2;
            else
                nparam = 1;
            end
            
            if(nparam ==1)
                X = log(Elem2param)/2;
            elseif(nparam == 2)
                X = zeros(nparam * nelem);
                X(1:2:end) = log(Elem2param(1,:))/2;
                X(2:2:end) = log(Elem2param(2,:))/2;
            else
                warning('the map from Elem2PAram to X is not unique!')
                X = zeros(nparam * nelem);
                X(1:3:end) = log(Elem2param(1,1,:))/2;
                X(3:3:end) = Elem2param(1,2,:) ./ sqrt(Elem2param(1,1,:));
                X(2:3:end) = log(Elem2param(2,2,:))/2 + log(Elem2param(1,1,:))/2 - log(Elem2param(1,2,:))/2;
            end
        end
        
        function lambda = getlformXe(Xe)
            nparam = length(Xe);
            if(nparam == 1)
                lambda = exp(2*Xe(1));
            elseif(nparam == 2)
                lambda = zeros(2);
                lambda(1,1) = exp(2*Xe(1));
                lambda(2,2) = exp(2*Xe(2));
            else
                lambda = zeros(2);
                lambda(1,1) = exp(2*Xe(1));
                lambda(2,2) = exp(2*Xe(2)) + Xe(3)^2;
                lambda(1,2) = exp(X(1)) * X(3);
                lambda(2,1) = lambda(1,2);
            end
        end
        
        
        function e2p = getElem2ParamfromX(X,nparam)
            nelem = length(X)/nparam;
            if(nparam == 1)
                e2p = exp(2*X);
            elseif(nparam == 2)
                e2p = zeros(2,nelem);
                e2p(1,:) = exp(2*X(1:2:end));
                e2p(2,:) = exp(2*X(2:2:end));
            else
                e2p = zeros(2,2,nelem);
                e2p(1,1,:) = exp(2*X(1:3:end));
                e2p(2,2,:) = exp(2*X(2:3:end)) + X(3:3:end)^2;
                e2p(1,2,:) = exp(X(1:3:end)) * X(3:3:end);
                e2p(2,1,:) = e2p(1,2,:);
            end
        end
        
    end
    
end

