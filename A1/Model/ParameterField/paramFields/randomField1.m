classdef randomField1 < pdeParameters
    %DETERMINITICFIELD1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % maindir is a matrix with as many columns as main directions (directions
        % of variation) of the random field and as many rows as field dimensions.
        maindir;
        
        % mainls is a vector of lengthscales corresponding to the main directions
        % given in maindir. It has as many entries as there are main directions.
        mainls;
        
        % Prefactor of covariance kernel
        sigma2;
        
        Info;
    end
    
    methods
        % Constructor
        function this = deterministicField1(lambdas, cps, rhos, maindir, mainls, sigma2)
            this = this@pdeParameters(lambdas, cps, rhos);
            
            this.Maindir = maindir;
            this.Mainls = mainls;
            this.Sigma2 = sigma2;
            
            this.Info = struct;
            this.Info.Type = class(this);
            this.Info.Maindir = maindir;
            this.Info.Mainls = mainls;
            this.Info.Sigma2 = sigma2;
        end
        
        function [ lambdas ] = evaluateParameterField(this, x_coords,y_coords)
            % maindir is a matrix with as many columns as main directions (directions
            % of variation) of the random field and as many rows as field dimensions.
            %
            % mainls is a vector of lengthscales corresponding to the main directions
            % given in maindir. It has as many entries as there are main directions.
            
            [row,col] = size(x_coords);
            
            if(row > 1)
                error('coordinate inputs to deterministicField1 have wrong format')
            end
            
            %Generate approximate Gaussian process sample functions in analytical form using Bochner's theorem
            nBasisFunctions = 1e3;
            
            % normalize columns of maindir and apply corresponding lengthscales
            for c = 1:length(this.Maindir(1,:))
                Lambda = this.Maindir(:,c)/norm(this.Maindir(:,c))/this.Mainls(c);
            end
            
            M = Lambda * Lambda';
            
            %Stacked samples from W
            W = mvnrnd(zeros(1, 2), M, nBasisFunctions);
            
            %Stacked samples from b
            b = repmat(2*pi*rand(nBasisFunctions, 1),[1,nx*ny]);
            
            %Draw coefficients gamma
            gamma = normrnd(0, 1, 1, nBasisFunctions);
            
            %Handle to sample function
            sampleFun = @(x) sqrt((2*this.Sigma2)/nBasisFunctions)*(gamma*cos(W*x + b));
            
            scal_lambdas = sampleFun([x_coords;y_coords]);
            scal_lambdas(scal_lambdas > 0) = this.Lambdas{1};
            scal_lambdas(scal_lambdas <= 0) = this.Lambdas{2};
            
            
            lambdas = zeros(2,2,col);
            lambdas(1,1,:) = scal_lambdas;
            lambdas(2,2,:) = scal_lambdas;
            
        end
        
    end
    
end

