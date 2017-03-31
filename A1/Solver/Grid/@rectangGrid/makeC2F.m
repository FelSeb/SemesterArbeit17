 function c2f = makeC2F(this, grid2)
            % Find out which grid is the coarse and which is the fine one
            if(this.N_elem > grid2.N_elem)
                cg = grid2;
                fg = this;
            elseif(this.N_elem < grid2.N_elem)
                fg = grid2;
                cg = this;
            else
                error('Both grids have the same number of elements. Coarse graining not possible.')
            end
            
            % check whether the Grids have the same geometry
            if( ~(isequal(cg.Geo_node{1}([1,end],[1,end]), fg.Geo_node{1}([1,end],[1,end])) && ...
                    isequal(cg.Geo_node{2}([1,end],[1,end]), fg.Geo_node{2}([1,end],[1,end]))) )
                error('The two grids do not have the same geometry. Coarse graining not possible.')
            end
            
            for ce = 1:cg.N_elem
                c2f{ce} = [];
                polygon = cg.Node2coord(:,cg.Elem2node(:,ce));
                % identify fg elements belonging to current coarse element
                for fe1 = 1:length(fg.Topo_elem(:,1))
                    for fe2 = 1:length(fg.Topo_elem(1,:))
                        
                        
                        
                        fe_x = fg.Geo_elem{1}(fe1,fe2);
                        fe_y = fg.Geo_elem{2}(fe1,fe2);
                        
                        isinpol = inpolygon(fe_x, fe_y,polygon(1,:),polygon(2,:));
                        
                        if(isinpol)
                            c2f{ce} = [c2f{ce}; fg.Topo_elem(fe2,fe1)];
                        end
                    end
                end
            end
        end