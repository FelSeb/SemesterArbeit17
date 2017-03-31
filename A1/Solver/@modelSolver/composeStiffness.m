function [this] = composeStiffness(this,varargin)
% Usage:
% [this] = composeStiffness(this)

E2Pdim = size(this.Grid.Elem2param);
if(length(E2Pdim) == 3)
    nparam = 3;
elseif(E2Pdim(1) == 2)
    nparam = 2;
else
    nparam = 1;
end

if(numel(varargin)>0)
    compute_grad = varargin{1};
else
    compute_grad = 0;
end

npe = this.Grid.FE.Nodes_per_elem;
nel = this.Grid.N_elem;
nno = this.Grid.N_node;

% Initialize finite element matricies
elem_mats = zeros(npe,npe,nel);
this.C = sparse(nno,nno);

% In case the derivatives of the finite element matix with respect to the
% material parameters are to be computed aswell...
if(compute_grad)
    elem_mats_grad = zeros(npe,npe,nel, nparam);
    this.C_grad = cell(1,nparam*nel);
    this.C_grad(:) = {sparse(nno,nno)};
else
    % Apparently this is needed so the pafor loop would work.
    elem_mats_grad = zeros(npe,npe,nel, nparam);
end

% parfor possible: This slows down the computation of the element matricies
% for small system sizes! Use parfor only for data generation!
for elem  = 1:nel
    % the coordinates of the element nodes
    xe = this.Grid.Node2coord(:,this.Grid.Elem2node(:,elem));
    
    % get local material parameter (distiguish 3 different formats)
    if(nparam == 3)
        lambda = this.Grid.Elem2param(:,:,elem);
    elseif(nparam == 2)
        lambda = diag(this.Grid.Elem2param(:,elem));
    else
        lambda = this.Grid.Elem2param(elem);
    end
    
    % integral (Gauss quadrature)
    for pxi = 1:length(this.Grid.FE.NumInt.in{1})
        for peta = 1:length(this.Grid.FE.NumInt.in{2})
            J = xe * squeeze(this.DN_all(pxi,peta,:,:))';
            
            detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1);
            invJ = 1/detJ * [J(2,2), -J(1,2); -J(2,1), J(1,1)];
            
            % compute coordinates
%             x_loc = xe * squeeze(this.N_all(pxi,peta,:));
%             x = x_loc(1);
%             y = x_loc(2);
            
            % Compose element matrix
            DN_A = squeeze(this.DN_all(pxi,peta,:,:));
            DN_B = squeeze(this.DN_all(pxi,peta,:,:));
            
            invJtDN_A = (invJ' * DN_A);
            invJtDN_B = (invJ' * DN_B);
            % integrand = (invJ' * DN_A)' * lambda * (invJ' * DN_B)
            % 'Diffusive matrix'
            elem_mats(:,:,elem) = elem_mats(:,:,elem) + ...
                this.Grid.FE.NumInt.weights(pxi,peta) * invJtDN_A' * lambda * invJtDN_B * detJ; % Cartesian
            
            % Compute element derivative matircies
            if(compute_grad)
                inds = (elem-1)*nparam+1 : (elem-1)*nparam+nparam;
                Xe = this.Grid.X(inds);
                elem_mats_grad(:,:,elem,:) = elem_mats_grad(:,:,elem,:) + ...
                    this.Grid.FE.NumInt.weights(pxi,peta) * detJ * ...
                    getElemMatGrad(Xe, invJtDN_A, invJtDN_B, this.Param_trafo);
            end
            
        end
    end
end

% Compose global matrix
for elem  = 1:nel
    % determine indices of matrix entries
    lins = this.Grid.Node2eq(this.Grid.Elem2node(:,elem));
    cols = this.Grid.Node2eq(this.Grid.Elem2node(:,elem));
    
    % 'Diffusive matrix'
    this.C(lins,cols) = this.C(lins,cols) + elem_mats(:,:,elem);
    
    % Compute derivative matircies
    if(compute_grad)
        inds = (elem-1)*nparam+1 : (elem-1)*nparam+nparam;
        for i = 1:length(inds)
            % Remeber: for nparam = 3, X has the following structure:
            %  X = [x1_11, x1_22, x1_12, x2_11, ...]'
            this.C_grad{inds(i)}(lins,cols) = elem_mats_grad(:,:,elem,i); 
        end
    end
    
end



end

function dfdX = getElemMatGrad(Xe, invJtDN_A, invJtDN_B, paramTrafo)

nparam = length(Xe);
if(nparam == 1)
    dfdl{1} = invJtDN_A(1,:)' * invJtDN_B(1,:) + invJtDN_A(2,:)' * invJtDN_B(2,:);
elseif(nparam == 2)
    dfdl{1} = invJtDN_A(1,:)' * invJtDN_B(1,:);
    dfdl{2} = invJtDN_A(2,:)' * invJtDN_B(2,:);
else
    dfdl{1} = invJtDN_A(1,:)' * invJtDN_B(1,:);
    dfdl{2} = invJtDN_A(2,:)' * invJtDN_B(2,:);
    dfdl{3} = (invJtDN_A(1,:)' * invJtDN_B(2,:) + invJtDN_A(2,:)' * invJtDN_B(1,:));
end

dldX = paramTrafo.getDldXe(Xe);

dfdX = zeros(4,4,nparam);
for ix = 1:nparam
    for il = 1:nparam
        dfdX(:,:,ix) = dfdX(:,:,ix) + dfdl{il} * dldX(il,ix);
    end
end

end




