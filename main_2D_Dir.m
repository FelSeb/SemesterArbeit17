clear;
%clc;
close all;

addpath('Util')
addpath('Anisotropic')
% This file produces a numerical 2D solution for the non-stationarry heat
% equation using a finite element approach.


% for convergence analysis:
nmethods = 1;
NNODES = [51,31,201,401];

for method_ind = 1:nmethods
    clearvars -except nmethods DT NNODES method_ind nmethods t05
    
    % Define parameters of the model
    s = SetModel();
    
    % Define parameters of the method
    nnodesx = NNODES(method_ind);
    nnodesy = NNODES(method_ind);
    m = SetMethod(s, nnodesx, nnodesy);
    
    % homogenize parameters
    [xmesh,ymesh] = meshgrid(linspace(0,s.domain(1,2),256),linspace(0,s.domain(2,2),256));
    [s.kappa_mat, s.kappa_mean] = Homogenize(xmesh,ymesh, s);

    % Compose global FE-matricies
    tic;
    C = ComposeC(s,m,nnodesx,nnodesy );
    toc;
    fprintf('Matrix composition complete.\n');
    
    % Apply Dirichlet BCs
    C = C(m.innereq,:);
    rhs = zeros(length(m.innereq),1) - C(:,m.boundeq_dir_all) * m.dir_val;
    C = C(:,m.innereq);
    
    % Solve linear system of equations
    tic;
    YOUT = C\rhs;
    toc;
    fprintf('System solved.\n');
    
    % post processing:
    [T, DTr, DTz] = PostProcessing( YOUT, nnodesx, nnodesy, m, method_ind);
    
    % movie-like plot of result
    SolPlot( T, DTr, DTz, s, m )
    
    % Interpret solution
    sol = MakeSolStruct(T, DTr, DTz, s, m);
    
end
%==========================================================================


% TO DO:
% -Try approch with power instead of intensity spectrum (Not a lot of benefit)
% -incorporate frequency/intesity/power info in the choice of lamnda_low and
% lambda_high
% - apply and test in problem
% - make matrix composition faster DONE
% - transform homogenized lambda tensor DONE
% - compare homogenized and non-homogenized results
% -compare coarse-grained homogenized (lamda tensor), coarse-grained
% homogenized (lambda scalar), fine-grained homogenized (lamda tensor),
% fine-grained homogenized (lambda-scalar), fine-grained heterogeneous
% -Problem Heat fluxes over all sides do not sum up to 0! -> only slow
% convergence: Testting showed that the total heatflux is roughly halfed if
% the number of nodes is quadrupeled
% -Modify code so that function for material parameter field can be
% exchanged easily...
% -Test matrix: 1.) Use heterogeneous field 2.) Use homogenized field 
% (lamda scalar mean) 3.) Use homogenized field (lambda tensor)
% -Display the difference between the solution using a heterogeneous
% parameter field and the one using a mean parameter field. If the
% difference is not significant the problem as such is not really
% interesting. The problem formulation would then need to be changed so
% that significant changes/ improvements would become visible.


% MOST IMPORTANT: Implement Framework!!!




