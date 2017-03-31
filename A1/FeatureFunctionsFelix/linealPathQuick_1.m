function [ features ] = linealPathQuick_1(lambda, param_trafo, direction, Grid, ce_node_coords, local_elem_numbers, opts)
% Compute the lineal path function for a linesegment of length z.
% ce_node_coords is a 2 x 4 matrix of the four nodes of the coarse element

if(direction > pi)
    error('linealPath function only allows directions between 0 and pi')
end

if(~isa(Grid, 'rectangGrid'))
    error('linealPathQuick_1 only defined for grid objects of type rectangGrid')
end

% low and high conductivity
lambda_hi = max(lambda);
lambda_lo = min(lambda);

if strcmp(opts.matrix_phase, 'lo')
    phase1 = lambda_lo;
elseif strcmp(opts.matrix_phase, 'hi')
    phase1 = lambda_hi;
else
    error('linealPathQuick_1 for high or low conducting phase as phase1?')
end

[sub_maps] = Grid.getSubMaps(local_elem_numbers);
sub_maps.Elem2param = lambda;

bins = zeros(1,opts.nbins);

% nlines = 100;

l = 1;
absl = 1;
while(l <= opts.nlines && absl <= opts.ntries)
    % choose as support points for the oriented lines multiple random convex
    % combinations of the ce_node_ccords
    r = 0.001+rand(4,1)*0.999;
    r = r./sum(r);
    x = ce_node_coords * r;
    %x = [0.495,0.4955];
    % direction = pi/4;
    
    
    segmentlength = getSegmentLengthQuick(x,direction,sub_maps,phase1);
    
    % Continue if randomly chosen poin does not lie in phase 1
    if(segmentlength == 0)
        absl = absl + 1;
        continue;
    end
    
    if(segmentlength < opts.z_max)
        ind = ceil(segmentlength/opts.z_max*opts.nbins);
    else
        ind = opts.nbins;
    end
    
    bins(1:ind)  = bins(1:ind) + 1 ;
    l = l+1;
    absl = absl + 1;
    %fprintf('segment finished \n')
end

% normalize histogram:
if(sum(bins)>0)
    bins = bins / sum(bins);
end

% Depending on specifications in opts transfrom output or not
if(strcmp(opts.propTo,'X'))
    features = bins;
elseif(strcmp(opts.propTo,'lambda'))
    features = zeros(opts.nbins,1);
    for ib = 1:opts.nbins
        features(ib) = param_trafo(bins(ib));
    end
else
   error('Missing specification in linealPathQuick_1 regarding parameter transformation') 
end



end


function [segmentlength] = getSegmentLengthQuick(x,dir,Grid, phase1)

%1.) Find element to which x belongs

%1.1) define a subset of all elements to which x may belong

%1.1.1) find the element of which the center is closest to x
dists = sqrt((Grid.Geo_elem{1}-x(1)).^2 + (Grid.Geo_elem{2}-x(2)).^2);
[min_dist, min_ind] = min(dists);
[min_min_dist, min_min_ind] = min(min_dist);
row = min_ind(min_min_ind);
col = min_min_ind;

neighbour = Grid.Topo_elem(row, col);

%2.) iterate over the neighbours that are determined by the line and
% compute the segment length from the distance between the initial
% element and the last element belongig to phase 1.
segmentlength = 0;
NDx = 0;
NDy = 0;
init_elem_coord = [Grid.Geo_elem{1}(row,col),Grid.Geo_elem{2}(row,col)];
while(phase1 == Grid.Elem2param(neighbour))
    final_elem_coord = [Grid.Geo_elem{1}(row,col),Grid.Geo_elem{2}(row,col)];
    
    %2.1) Get next element
    [neighbour, row, col, NDx,NDy] = getNeighbour(row, col, NDx, NDy, dir, Grid);
    
    % Stop if domain is exceeded
    if(neighbour ==-1 || Grid.Elem2param(neighbour) ~= phase1)
        segmentlength = norm(init_elem_coord-final_elem_coord);
        break;
    end
end


end


function [neighbour, row, col, NDx_,NDy_] = getNeighbour(row, col, NDx,NDy, dir, Grid)
% In this function it is assumed that all elements have the same size and
% shape. It should only be applied if fine element size is sufficiently
% small.

if(dir < pi/2)
    dir_approx(1) = (NDx+1)/sqrt((NDx+1)^2+NDy^2);
    dir_approx(2) = (NDx+1)/sqrt((NDx+1)^2+(NDy+1)^2);
else
    dir_approx(1) = (NDx-1)/sqrt((NDx-1)^2+(NDy)^2);
    dir_approx(2) = (NDx-1)/sqrt((NDx-1)^2+(NDy+1)^2);
end
dir_approx(3) = NDx/sqrt(NDx^2+(NDy+1)^2);


delta_angle = abs(dir-acos(dir_approx));
[min_angle, ind_min] = min(delta_angle);

if(ind_min == 3)
    NDy_ = NDy+1;
    NDx_ = NDx;
elseif(ind_min == 2)
    if(dir < pi/2)
        NDx_ = NDx+1;
        NDy_ = NDy+1;
    else
        NDx_ = NDx-1;
        NDy_ = NDy+1;
    end
else
    if(dir < pi/2)
        NDx_ = NDx+1;
        NDy_ = NDy;
    else
        NDx_ = NDx-1;
        NDy_ = NDy;
    end
end

col = col + (NDx_-NDx);
row = row + (NDy_-NDy);

if(row > Grid.Ny ||row < 1 || col > Grid.Nx || col < 1)
    neighbour = -1;
else
    neighbour = Grid.Topo_elem(row,col);
end


end





