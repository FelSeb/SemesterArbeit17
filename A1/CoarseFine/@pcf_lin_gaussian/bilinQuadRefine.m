function [ fine_field ] = bilinQuadRefine(this, coarse_field )

% Bilinear interpolation 

% Values of the four corresponding nodes
vals = [coarse_field(1,1);
    coarse_field(1,2);
    coarse_field(2,2);
    coarse_field(2,1);];

interpolated_value = this.Offset + this.Interpolation_matrix * vals;

fine_field = reshape(interpolated_value,size(this.Fine_grid{1}));

figure;
surf(this.Fine_grid{1},this.Fine_grid{2},fine_field)
view(0,90);

end

