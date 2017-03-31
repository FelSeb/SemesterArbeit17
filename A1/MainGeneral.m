clear;
clc;
close all;

addpath('Model')
addpath('Solver')
addpath('Solver/Element')
addpath('Solver/Grid')
addpath('Model/ParameterField')
addpath('Model/Domain')

My_pdeparams = deterministicField1({1}, {1}, {1});

My_domain = rectangDomain(0,1,0,1);

My_bound_conds = boundaryConds(My_domain);
My_bound_conds = My_bound_conds.SetUpPolynimialBoundConds([0,0,0,1]);

My_PrimaryModel = primaryModel(My_domain, My_pdeparams, My_bound_conds);

My_NumericModel = modelSolver(My_PrimaryModel);

My_Grid = rectangGrid(2,2,My_PrimaryModel,'quadrilateral','sf_lin','GaussQuad2');

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

My_NumericModel.plotSol(0);












