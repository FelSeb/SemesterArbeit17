clear all;
close all;
%clc;

addpath('CoarseFine')
addpath('CoarseGrain')
addpath('Model')
addpath('Solver')
addpath('Solver/Element')
addpath('Solver/Grid')
addpath('Model/ParameterField')
addpath('Model/ParameterField/paramFields')
addpath('Model/Domain')
addpath('Training')
addpath('Training/MCMC')
addpath('Training/ProposalDensities')
addpath('featureTrafo')
addpath('FeatureFunctionsFelix')
addpath('Training/VariationalBayes')

warning('use genpath funciton here!!!')
%% Set parallel computing options
run_parallel = 1;
if(numel(gcp('nocreate')) ==0 && run_parallel)
    nworkers = 4;
    parpool(4);
elseif(numel(gcp('nocreate')) ~= 0 && ~run_parallel)
    poolobj = gcp('nocreate');
    delete(poolobj);
end

%% Primary Model
% Remember: anisotropic random fields are not considered!
%My_pdeparams = deterministicField1_scalar({1}, {1}, {1});
maindir1 = [1,0];
maindir2 = [0,1];
volume_fraction = 0.35; % V_hi/V_ges;
My_pdeparams = randomField1_scalar({1,10}, {1,10}, {1,10}, volume_fraction,[maindir1',maindir2'], [0.08,0.08], 1);

% Plot an exemplary realization of the parameter field
[xg,yg] = meshgrid(linspace(0,1,50),linspace(0,1,50));
% My_pdeparams.plotField(xg,yg);

My_domain = rectangDomain(0,1,0,1);
My_bound_conds = boundaryConds(My_domain);
My_bound_conds = My_bound_conds.SetUpPolynimialBoundConds([0,150,100,-30]);

My_PrimaryModel = primaryModel(My_domain, My_pdeparams, My_bound_conds);

%% Coarse and fine grids
Coarse_grid = rectangGrid(2,2,My_PrimaryModel,'quadrilateral','sf_lin','GaussQuad2' );
Fine_grid = rectangGrid(32,32,My_PrimaryModel,'quadrilateral','sf_lin','GaussQuad2' );
Coarse_grid = Coarse_grid.makeMaps(My_PrimaryModel);
Fine_grid = Fine_grid.makeMaps(My_PrimaryModel);

c2f = Coarse_grid.makeC2F(Fine_grid);

%% Initialize coarse and fine solvers
param_trafo = paramTrafoPoDe();
coarseSolver = modelSolver(My_PrimaryModel, param_trafo);
fineSolver = modelSolver(My_PrimaryModel);

%% PCF
init_var = 1e5; % This parameter must be set to relatively high value (maybe in relation to the temperature values that occur)
My_pcf = pcf_lin_gaussian(Fine_grid,Coarse_grid,c2f,'varS', coarseSolver, 1, init_var);
My_pcf = My_pcf.initParameters;

%% PC
init_var = 1;
FeatFun = @(lambdas) featfun3_1(lambdas,c2f, param_trafo, Fine_grid, Coarse_grid); % make feature function a seperate class
My_pc = pc_lin_feature_gaussian1(Fine_grid,Coarse_grid,c2f,FeatFun,init_var);
My_pc = My_pc.initParameters;

%% Data Generator
n_data_points = 30; 
tic;
generator = DataGenerator(Fine_grid, fineSolver, 'Data1/Dataset32');
N_batch = 1;
generator.generate(n_data_points,N_batch);
fprintf('Data generation complete\n');
toc;

%% Data Provider
provider = DataProvider('Data1/Dataset32');

%% For test purposes display generated data
% for i = 1:5
%     [x_i, y_i] = provider.provideDataPoint(i);
%     PostSol.T = y_i;
%     Fine_grid.plotSol(PostSol, 0);
% end

%% mhMCMC
prop_var = 0.6;
My_proposer = MALAProp(prop_var);
% prop_var = 5e-1;
% My_proposer = GaussianProp(prop_var);

n_MC_samples = 33;
My_Sampler = mhMCMC(My_proposer, n_MC_samples);

%% EM
n_training_samples = n_data_points;
maxiter = 5;
My_EM = EM(My_pc, My_pcf, My_PrimaryModel, n_data_points, provider, maxiter, My_Sampler);
[theta_pc, theta_pcf, My_EM] = My_EM.optimizeEM;


%% Display the fitted values 


%%
% To Do / Questions
% 1.) Parallelism should be controlable ...remote control whether parfor is
% used or not.

% 2.) Implement more feature function. Also for the anisotropic case...

