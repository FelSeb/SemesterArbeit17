clear all;
close all;
clc;

addpath('CoarseFine')
addpath('Model')
addpath('Solver')
addpath('Solver/Element')


My_pdeparams = pdeParameters({1}, {1}, {1}, 'deterministic1');

My_domain = domain();
My_domain = My_domain.rectangularDomain(0,1,0,1);

My_bound_conds = boundaryConds(My_domain);
My_bound_conds = My_bound_conds.SetUpPolynimialBoundConds([0,0,0,1]);

My_PrimaryModel = primaryModel(My_domain, My_pdeparams, My_bound_conds);


Coarse_grid = rectangGrid(2,2,My_PrimaryModel,'quadrilateral','sf_lin','GaussQuad2' );
Fine_grid = rectangGrid(4,4,My_PrimaryModel,'quadrilateral','sf_lin','GaussQuad2' );
Coarse_grid = Coarse_grid.makeMaps(My_PrimaryModel);
Fine_grid = Fine_grid.makeMaps(My_PrimaryModel);

c2f = Coarse_grid.makeC2F(Fine_grid);

My_pcf = pcf_lin_gaussian(Fine_grid,Coarse_grid,c2f,'allfix', 1);

My_pcf = My_pcf.getInterpMat(1);


% Test interpolation matrix:
coarse_vals = My_pcf.Coarse_grid.Node2coord(1,1:My_pcf.Coarse_grid.N_node)' .* ...
    My_pcf.Coarse_grid.Node2coord(2,1:My_pcf.Coarse_grid.N_node)';

figure;
surf(My_pcf.Coarse_grid.Geo_node{1}, My_pcf.Coarse_grid.Geo_node{2}, coarse_vals(My_pcf.Coarse_grid.Topo_node));

fine_vals = My_pcf.Interpolation_matrix * coarse_vals;
figure;
surf(My_pcf.Fine_grid.Geo_node{1}, My_pcf.Fine_grid.Geo_node{2}, fine_vals(My_pcf.Fine_grid.Topo_node));















