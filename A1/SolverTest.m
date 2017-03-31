

clear;
clc;
close all;

addpath('Model')
addpath('Solver')
addpath('Solver/Element')
addpath('Solver/Grid')
addpath('Model/ParameterField')
addpath('Model/Domain')

% Using a polynomial function as analytical reference solution:

% Array of coefficients:
coeffs = [1,1,1,1];

%% Compute numerical solution
My_pdeparams = deterministicField1({1}, {1}, {1});
% maindir1 = [1,0];
% maindir2 = [0,1];
% volume_fraction = 0.35; % V_hi/V_ges;
% My_pdeparams = randomField1_scalar({1,10}, {1,10}, {1,10}, volume_fraction,[maindir1',maindir2'], [0.08,0.08], 1);


My_domain = rectangDomain(0,1,0,1);

My_bound_conds = boundaryConds(My_domain);
My_bound_conds = My_bound_conds.SetUpPolynimialBoundConds(coeffs);

My_PrimaryModel = primaryModel(My_domain, My_pdeparams, My_bound_conds);

My_NumericModel = modelSolver(My_PrimaryModel);

My_Grid = rectangGrid(32,32,My_PrimaryModel,'quadrilateral','sf_lin','GaussQuad2');

My_NumericModel = My_NumericModel.setGrid(My_Grid);

My_NumericModel = My_NumericModel.preEval;

tic;
My_NumericModel = My_NumericModel.composeStiffness;
toc;

My_NumericModel = My_NumericModel.applyBoundCond;

tic;
My_NumericModel = My_NumericModel.solveSystem;
toc;

tic;
My_NumericModel = My_NumericModel.getPostSol(0);
toc;

%% Evaluate analytical reference solution:
ref_sol = My_bound_conds.evaluateBoundaryFun(My_NumericModel.Grid.Geo_node{1},My_NumericModel.Grid.Geo_node{2},1);

%% Plot 
My_NumericModel.plotSol(0);
title('Computed')
figure;
surf(My_NumericModel.Grid.Geo_node{1},My_NumericModel.Grid.Geo_node{2},ref_sol);
title('Reference')

%% Compute error:
formated_ref_sol = reshape(ref_sol',size(My_NumericModel.PostSol.T));
error = formated_ref_sol - My_NumericModel.PostSol.T;

PostSol.T = error;
My_Grid.plotSol(PostSol, 0)
xlabel('x')
ylabel('y')
zlabel('Delta T')
title('Error')

if(max(max(abs(error))) < 1e-10)
   fprintf('Unit test succesfull') 
end






