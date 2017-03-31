function [ features ] = featfun3_1( lambdas , c2f, param_trafo, Fine_grid, Coarse_grid)
% This functions returns the feature matrix. Each column contains one
% feature in vector format. 
% lambdas       vector of scalar material parameters for all fine elements
% c2f           coarse elements to fine elements map
% 
% features      feature matrix. Each column is one feature. One entry of a
%               column of the feature matrix is associated with one
%               coarse element.

% lambdas has Elem2param format or vector format. Since the material 
% parameters of the fine field are scalar, the two formats are identical.

N_features = 1;
N_coarse_elem = numel(c2f);

features = zeros(N_coarse_elem, N_features);


for ce = 1:N_coarse_elem
    % Feature 1: local mean
    lambda_eff = mean(lambdas(c2f{ce}));
    features(ce,1) = param_trafo.getXfromElem2param(lambda_eff);
    
    % Feature 2.1: Maxwell-Garnett, high cond matrix
    features(ce,2) = maxwellGarnett_1(lambdas(c2f{ce}), param_trafo , 'hi');
    
    % Feature 2.2: Maxwell-Garnett, lo cond matrix
    features(ce,3) = maxwellGarnett_1(lambdas(c2f{ce}), param_trafo , 'lo');
    
    % Feature 3.1: DEM, high cond matrix
    features(ce,4) = DEM_1(lambdas(c2f{ce}), param_trafo , 'hi');
    
    % Feature 3.2: DEM, lo cond matrix
    features(ce,5) = DEM_1(lambdas(c2f{ce}), param_trafo , 'lo');
    
    % Feature 3.2: DEM, lo cond matrix
    features(ce,6) = SCA_1(lambdas(c2f{ce}), param_trafo);
    
%     % Feature 4: liealPath
%     opts.nbins = 5;
%     opts.nlines = 50;
%     opts.z_max = 0.5;
%     opts.propTo = 'X';
%     opts.ntries = 2;
%     n_linpath_feat = opts.nbins;
%     ce_node_coords = Coarse_grid.Node2coord(:,Coarse_grid.Elem2node(:,ce));
%     
%     % Feature 4.1: liealPath, angle: 0
%     opts.matrix_phase = 'hi';
%     features(ce,6:6+n_linpath_feat-1) = linealPathQuick_1(lambdas(c2f{ce}), param_trafo, 0, Fine_grid, ce_node_coords, c2f{ce}, opts);
% 
%     % Feature 4.2: liealPath, angle: pi/2
%     opts.matrix_phase = 'hi';
%     features(ce,6+n_linpath_feat:6+2*n_linpath_feat-1) = linealPathQuick_1(lambdas(c2f{ce}), param_trafo, pi/2, Fine_grid, ce_node_coords, c2f{ce}, opts);
% 
%     % Feature 4.3: liealPath, angle: 0
%     opts.matrix_phase = 'lo';
%     features(ce,6+2*n_linpath_feat:6+3*n_linpath_feat-1) = linealPathQuick_1(lambdas(c2f{ce}), param_trafo, 0, Fine_grid, ce_node_coords, c2f{ce}, opts);
% 
%     % Feature 4.4: liealPath, angle: pi/2
%     opts.matrix_phase = 'lo';
%     features(ce,6+3*n_linpath_feat:6+4*n_linpath_feat-1) = linealPathQuick_1(lambdas(c2f{ce}), param_trafo, pi/2, Fine_grid, ce_node_coords, c2f{ce}, opts);

end


end






