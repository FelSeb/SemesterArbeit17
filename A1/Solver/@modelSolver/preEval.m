function [ this ] = preEval( this )

% Some precomputation to improve performance
nxi = length(this.Grid.FE.NumInt.in{1}); 
neta = length(this.Grid.FE.NumInt.in{2});
this.N_all = zeros(nxi,neta,this.Grid.FE.Nodes_per_elem);
this.DN_all = zeros(nxi,neta,2,this.Grid.FE.Nodes_per_elem);

for  node = 1:this.Grid.FE.Nodes_per_elem % Shape
    for pxi = 1:nxi
        for peta = 1:neta
            [this.N_all(pxi,peta,node),this.DN_all(pxi,peta,:,node)] = ...
                this.Grid.FE.SF(this.Grid.FE.NumInt.in{1}(pxi,peta),this.Grid.FE.NumInt.in{2}(pxi,peta), node);
        end
    end
end

% Some precomputation to improve performance (1D)
nxi1D = length(this.Grid.FE.NumInt1D.in);
this.N_all1D = zeros(nxi,this.Grid.FE.Nodes_per_elem1D);
this.DN_all1D = zeros(nxi,this.Grid.FE.Nodes_per_elem1D);

for  node = 1:this.Grid.FE.Nodes_per_elem1D % Shape
    for pxi = 1:nxi1D
            [this.N_all1D(pxi,node),this.DN_all1D(pxi,node)] = ...
                this.Grid.FE.SF1D(this.Grid.FE.NumInt1D.in(pxi), node);
    end
end

end

