clear;
clc;
close all;

%this is a test.
%this is a second change.

addpath('Util')
% This file produces a numerical 1D solution for the non-stationarry heat
% equation using a finite element approach.
% We assume a constant boundary condition at x = L and a Neumann BC at
% x = 0.

% for convergence analysis:
nmethods = 1; % method 1 is already converged!
fi_mi = 2/10;
mi_co = 4/10;
NNODES_fi = [21,51,102];
NNODES_mi = [6,31,31];
NNODES_co = [6,31,31];

for method_ind = 1:nmethods
    clearvars -except nmethods DT ...
        NNODES_fi NNODES_mi NNODES_co fi_mi mi_co ...
        method_ind nmethods t05
    
    % initialize struct
    s = struct;
    
    % Define some parameters of the model
    s.rho = 100;
    s.cp = 100;
    s.lambda = 1.;%0.5;
    s.kap = s.lambda/(s.rho*s.cp);
    
    % define the heat influx (function) at x = 0;
    s.influx0 = 0.1;
    s.influx1 = -0.1;
    s.T0 = 5;
    
    % Define some parameters of the method
    domain = [0,50]; % Domain
    nodel = 2; % nodes per element (fixed)
    
    gps = [-1,1]/sqrt(3); % Gauss quadrature points
    gw = [1,1];% Gauss quadrature weights
    
    %     % Generate a non-uniform grid
    %     fi_co = 1/5;
    %     nnod_fi = NNODES_fi(method_ind);
    %     nnod_co = NNODES_co(method_ind);
    %     grid_fine = linspace(domain(1),domain(2)*fi_mi,nnod_fi);
    %     h_fine = grid_fine(2)-grid_fine(1);
    %     grid_middle = linspace(grid_fine(end) + h_fine ,domain(2)*mi_co,nnod_co);
    %     h_middle = grid_fine(2)-grid_fine(1);
    %     grid_coarse = linspace(grid_middle(end) + h_middle ,domain(2),nnod_co);
    %     grid = [grid_fine,grid_middle,grid_coarse];
    
    % Generate a uniform grid
    grid = linspace(domain(1),domain(2),11);
    
    nnodes = length(grid); % number of nodes
    nelems = nnodes - 1; % number of elements
    
    nod2coord = grid; % node to coordinate mapping
    elem2nod(1,:) = 1:nelems; % element to node mapping
    elem2nod(2,:) = 2:nelems+1; % element to node mapping
    nod2eq = zeros(1,nnodes);
    nod2eq(1:nnodes) = 1:nnodes;% node to equation mapping
    
    % Define shape- and test functions and their derivatives
    addpath('ShapeTest')
    
    % Initialize finite element matricies
    C = zeros(nnodes); % Diffusive matrix
    
    % Compose global matricies
    for elem  = 1:nelems
        
        % Compose element matrix
        for nodA = 1:nodel % Test
            for  nodB = 1:nodel % Shape
                
                % determine indices of matrix entries
                lin = nod2eq(elem2nod(nodA,elem));
                col = nod2eq(elem2nod(nodB,elem));
                
                % integral (Gauss quadrature)
                for gp = 1:length(gps)
                    h = nod2coord(elem2nod(2)) - nod2coord(elem2nod(1));
                    % jacobi determinant of transformation
                    J = 0.5 * h; % = dx/dxi (in general depending on xi)
                    invJ = 1/J; % = dxi/dx;
                    
                    % Evaluation of shape- and test funcitons
                    [N_A, DN_A] = shape_lin_1D(gps(gp), nodA);
                    [N_B, DN_B] = shape_lin_1D(gps(gp), nodB);
                    
                    % 'Diffusive matrix'
                    C(lin,col) = C(lin,col) + gw(gp) * s.kap * DN_A * invJ * DN_B * invJ * J;
                end
                
            end
        end
    end
    
%     % Neu Neu (doesn't make sense)
%     rhs = zeros(nnodes,1);
%     rhs(1) = s.kap/ s.lambda * s.influx0; % flux boundary condition at x = 0
%     rhs(end) = s.kap/ s.lambda * s.influx1;% temperature boundary condition at x = L
%     
%     YOUT = C\rhs;
%     
%     %plot
%     plot(grid,YOUT)
    
    % Neu Dir
    rhs = zeros(nnodes-1,1);
    rhs(1) = s.kap/ s.lambda * s.influx0; % flux boundary condition at x = 0
    rhs(end) = -C(end-1,end)*s.T0;% temperature boundary condition at x = L
    C = C(1:end-1,1:end-1);
    
    YOUT = C\rhs;
    
    %plot
    figure();
    plot(grid,[YOUT;s.T0])
     
end






















