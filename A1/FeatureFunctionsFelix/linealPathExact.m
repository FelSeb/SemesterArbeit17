function [ bins ] = linealPathExact(x_tensor, phase1, Grid, ce_node_coords, direction)
% Compute the lineal path function for a linesegment of length z.
% If x_tensor is of a higher order format each component is considered
% individually.

% ce_node_coords is a 2 x 4 matrix of the four nodes of the coarse element

if(direction > pi)
   error('linealPath function only allows directions between 0 and pi') 
end

% Set parameter field
Grid = Grid.setElem2param(x_tensor);

nbins = 10;
maxlen = 0.3;
bins = zeros(1,nbins);

nlines = 100;

l = 1;
while(l <= nlines)
    % choose as support points for the oriented lines multiple random convex
    % combinations of the ce_node_ccords
    r = 0.001+rand(4,1)*0.999;
    r = r./sum(r);
    x = ce_node_coords * r;
    %x = [0.495,0.4955];
    % direction = pi/4;
    
    
    segmentlength = getSegmentLength(x,direction,Grid,phase1);
    
    % Continue if randomly chosen poin does not lie in phase 1
    if(segmentlength == 0)
        continue;
    end
    
    if(segmentlength < maxlen)
        ind = ceil(segmentlength/maxlen*nbins);
    else
        ind = nbins;
    end
    
    bins(1:ind)  = bins(1:ind) + 1 ;
    l = l+1;
    fprintf('segment finished \n')
end

end



function angles = getElementAngles(elemnodecoords)

diagangle = atan( (elemnodecoords(2,3)-elemnodecoords(2,2))/...
    (elemnodecoords(1,2)-elemnodecoords(1,1)) );


angles = [diagangle, diagangle + pi-2*diagangle, pi + diagangle, 2*pi-diagangle];

end

function [intersectionsx,intersectionsy,neighbour1,neighbour2] = getElemSegment(x,dir,elemnodecoords)

x1 = [elemnodecoords(1,:),elemnodecoords(1,1)];
y1 = [elemnodecoords(2,:),elemnodecoords(2,1)];

slope = tan(dir);
linex = [min(x1)-1,max(x1)+1];
liney = x(2)+slope*(linex-x(1));

[intersectionsx,intersectionsy,ii] = polyxpoly(y1,x1,linex,liney);
neighbour1 = ii(1,1);
neighbour2 = ii(2,1);

end



function [segmentlength] = getSegmentLength(x,dir,Grid, phase1)
if(isa(Grid, 'rectangGrid'))
    %1.) Find element to which x belongs
    
    %1.1) define a subset of all elements to which x may belong
    
    %1.1.1) find the element of which the center is closest to x
    dists = sqrt((Grid.Geo_elem{1}-x(1)).^2 + (Grid.Geo_elem{2}-x(2)).^2);
    [min_dist, min_ind] = min(dists);
    [min_min_dist, min_min_ind] = min(min_dist);
    row = min_ind(min_min_ind);
    col = min_min_ind;
    
    rowup = min(row+1,Grid.Ny);
    rowlo = max(row-1,1);
    colup = min(col+1,Grid.Nx);
    collo = max(col-1,1);
    
    midelem = Grid.Topo_elem(row, col);
    
    %1.1.2) define subset
    elem_subset = [midelem,...
        Grid.Topo_elem(rowup, colup),...
        Grid.Topo_elem(rowup, col),...
        Grid.Topo_elem(rowup, collo),...
        Grid.Topo_elem(row, colup),...
        Grid.Topo_elem(row, collo),...
        Grid.Topo_elem(rowlo, colup),...
        Grid.Topo_elem(rowlo, col),...
        Grid.Topo_elem(rowlo, collo),...
        ];
    
    roco = [0,1,1,1,0,0,-1,-1,-1;
        0,1,0,-1,1,-1,1,0,-1];
    
    %1.2) find element to which x belongs
    for e = 1:9
        elemnodecoords = Grid.Node2coord(:,Grid.Elem2node(:,elem_subset(e)));
        x1 = [elemnodecoords(1,:),elemnodecoords(1,1)];
        y1 = [elemnodecoords(2,:),elemnodecoords(2,1)];
        if(inpolygon(x(1),x(2),y1,x1))
            xelem = elem_subset(e);
            row = row+roco(1,e);
            col = col+roco(2,e);
            break;
        end
    end
    
    %2.) iterate over the neighbours that are determined by the line and
    %add up the respective linesegments until a new phase is encountered.
    segmentlength = 0;
    while(phase1 == Grid.Elem2param(xelem))
        %2.1) Get linesegment and next element
        elemnodecoords = Grid.Node2coord(:,Grid.Elem2node(:,xelem));
        [intersectionsx,intersectionsy,n1,n2] = getElemSegment(x,dir,elemnodecoords);
        
        % make sure to move in one direction
        if((intersectionsx(1)>intersectionsx(2) && dir < pi/2 ||...
                intersectionsy(1)>intersectionsy(2) && dir >= pi/2))
            interx = [intersectionsx(1),intersectionsy(1)];
            neighbour = n1;
        else
            interx = [intersectionsx(2),intersectionsy(2)];
            neighbour = n2;
        end
        
        elem_segmentlength = sqrt((x(1)-interx(1)).^2 + (x(2)-interx(2)).^2);
        segmentlength = segmentlength + elem_segmentlength;
        
        if(neighbour == 1)
            col = col - 1;
        elseif(neighbour == 2)
            row = row + 1;
        elseif(neighbour == 3)
            col = col + 1;
        else
            row = row - 1;
        end
        
        % Stop if domain is exceeded
        if(col > size(Grid.Topo_elem,2) || row > size(Grid.Topo_elem,1))
            break;
        end
%         % Alternatively: periodic boundary conditions. There must be a
%         % max segment length then...
%         if(col > size(Grid.Topo_elem,2))
%             col = 1;
%         elseif(row > size(Grid.Topo_elem,1))
%             row = 1;
%         end
        
        xelem = Grid.Topo_elem(row,col);
        x = interx;
        
        nextphase = Grid.Elem2param(xelem);
        % Stop if phase 2 is encountered
        if(nextphase ~= phase1)
            break;
        end
    end
else
    error('function not defined for grids that are not of type rectangGrid')
end


end