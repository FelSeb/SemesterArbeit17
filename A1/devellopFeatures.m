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
addpath('Model/Domain')

addpath('MYFeatureFunctions')

nx = 256;
ny = 256;

%% Primary Model
%My_pdeparams = deterministicField1({1}, {1}, {1});
%My_pdeparams = deterministicField1_scalar({1}, {1}, {1});
maindir1 = [1,0];
maindir2 = [0,1];
maindirs = maindir1';%[maindir1',maindir2'];
lengthscales = 0.05;%[0.05,0.05];
My_pdeparams = randomField1_scalar({1,2}, {1,2}, {1,2}, maindirs, lengthscales, 1);

My_domain = rectangDomain(0,1,0,1);

My_bound_conds = boundaryConds(My_domain);
My_bound_conds = My_bound_conds.SetUpPolynimialBoundConds([0,0,0,1]);

My_PrimaryModel = primaryModel(My_domain, My_pdeparams, My_bound_conds);

%% Grid
Fine_grid = rectangGrid(nx,ny,My_PrimaryModel,'quadrilateral','sf_lin','GaussQuad2' );
Fine_grid = Fine_grid.makeMaps(My_PrimaryModel);


% Plot an exemplary realization of the parameter field and use it to test
% the linealPath function
params = My_pdeparams.plotField(Fine_grid.Geo_elem{1},Fine_grid.Geo_elem{2});
x_tensor = reshape(params,1,numel(params));


%% Feature test
% tic;
% bins = linealPathExact(x_tensor, 1, Fine_grid, [0,1,1,0; 0,0,1,1], atan(0.5))
% toc;


tic;
bins = linealPathQuick(x_tensor, 1, Fine_grid, [0,1,1,0; 0,0,1,1], pi/2)
toc;






