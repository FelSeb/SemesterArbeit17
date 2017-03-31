classdef BayesianCGPDE
    %BAYESIANCGPDE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        primaryModel; % Model that is to be analyzed
        
        modelSolver; % Numerical Algorithm to solve PDE
        
        pcf; % Probabilistic coarse2fine map
        
        pc; % Probabilitsic fine2coarse map
        
        n_training_samples; % Number of training samples
        
        path_to_training_data; % directory of training data
        
        path_to_MC_samples % directory of MC samples
    end
    
    methods
        
        
        % Genearte training data
        generateTrainingData(this);
        
        % Fit the parameters of pcf and pc
        fit(this);
        
        % Obrain predictions regarding the statistics of the output
        predict(this);
       
        
        % Check the functionality of the FEM by comparing it to an
        % analytical reference solution
        checkFEMsol(this)
        
        
        
        
        
        
    end
end

