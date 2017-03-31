clear all;
close all;
clc;

% This file illustrates how (anisotropic) Bochner samples may be generated

volume_fraction = 0.3; % V_hi/V_ges;
threshold = erfinv(volume_fraction *2 -1)*sqrt(2);

nx = 100;
ny = 100;

[xg,yg] = meshgrid(linspace(0,1,nx),linspace(0,1,ny));

coords = zeros(2,nx*ny);
coords(1,:) = reshape(xg,1,nx*ny);
coords(2,:) = reshape(yg,1,nx*ny);


% Anisotropic:


%Generate approximate Gaussian process sample functions in analytical form using Bochner's theorem
sigma_f2 = 1;
nBasisFunctions = 1e3;


dir_1 = [1,0];
dir_2 = [0,1];
%dir_3 = [1,1];
ls_1 = 0.05;
ls_2 = 0.05;
ls_3 = 0.01;

dim = 2;
n_main_dir = 2;
Lambda = zeros(dim,n_main_dir);
Lambda(:,1) = dir_1/norm(dir_1) / ls_1;
Lambda(:,2) = dir_2/norm(dir_2) / ls_2;
% Lambda(:,3) = dir_3/norm(dir_3) / ls_3;
M = Lambda * Lambda';

%Stacked samples from W, see reference_notes
W = mvnrnd(zeros(1, 2), M, nBasisFunctions);

%Stacked samples from b, see notes
b = repmat(2*pi*rand(nBasisFunctions, 1),[1,nx*ny]);


for i = 1:1
    %Draw coefficients gamma
    gamma = normrnd(0, 1, 1, nBasisFunctions);
    
    %Handle to sample function
    sampleFun = @(x) sqrt((2*sigma_f2)/nBasisFunctions)*(gamma*cos(W*x + b));

    res = sampleFun(coords);
    
%     var(reshape(res,1,numel(res)))
%     mean(reshape(res,1,numel(res)))
    
    res = reshape(res,nx,ny);
    res = (res < threshold)*1.0;
    
    vol_frac = sum(reshape(res,1,numel(res)))/numel(res)
    
    pcolor(xg,yg,res)
    
    %view([0,90])
    
    %pause();
end




