classdef randomField1_scalar < pdeParameters
    %DETERMINITICFIELD1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % maindir is a matrix with as many columns as main directions (directions
        % of variation) of the random field and as many rows as field dimensions.
        Maindir;
        
        % mainls is a vector of lengthscales corresponding to the main directions
        % given in maindir. It has as many entries as there are main directions.
        Mainls;
        
        % Prefactor of covariance kernel
        Sigma2;
        
        % Volume fraction of the high conducting material
        Volume_fraction; % V_hi/V_ges;
        
        Info;
    end
    
    methods
        % Constructor
        function this = randomField1_scalar(lambdas, cps, rhos, vol_frac, maindir, mainls, sigma2)
            this = this@pdeParameters(lambdas, cps, rhos);
            
            this.Volume_fraction = vol_frac;
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
            
            
            threshold = erfinv(this.Volume_fraction *2 -1)*sqrt(2);
            
            [row,col] = size(x_coords);
            
            if(row > 1)
                error('coordinate inputs to deterministicField1_scalar have wrong format')
            end
            
            %Generate approximate Gaussian process sample functions in analytical form using Bochner's theorem
            nBasisFunctions = 1e3;
            
            % normalize columns of maindir and apply corresponding lengthscales
            n_main_dir = length(this.Maindir(1,:));
            dim = 2;
            Lambda = zeros(dim,n_main_dir);
            for c = 1:length(this.Maindir(1,:))
                Lambda(:,c) = this.Maindir(:,c)/norm(this.Maindir(:,c))/this.Mainls(c);
            end
            
            M = Lambda * Lambda';
            
            %Stacked samples from W
            W = mvnrnd(zeros(1, 2), M, nBasisFunctions);
            
            %Stacked samples from b
            b = repmat(2*pi*rand(nBasisFunctions, 1),[1,col]);
            
            %Draw coefficients gamma
            gamma = normrnd(0, 1, 1, nBasisFunctions);
            
            %Handle to sample function
            sampleFun = @(x) sqrt((2*this.Sigma2)/nBasisFunctions)*(gamma*cos(W*x + b));
            
            
            if(this.Lambdas{1} > this.Lambdas{2})
               lambda_hi = this.Lambdas{1};
               lambda_lo = this.Lambdas{2};
            else
               lambda_hi = this.Lambdas{2};
               lambda_lo = this.Lambdas{1};
            end
            
            scal_lambdas = sampleFun([x_coords;y_coords]);
            lambdas = zeros(size(scal_lambdas));
            lambdas(scal_lambdas < threshold) = lambda_hi;
            lambdas(scal_lambdas >= threshold) = lambda_lo;
            
            
           
        end
        
    end
    
end

