classdef varDistGaussian
    %VARDIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Mean parameters
        Phi_mu;
        
        % Covariance parameters (must be transformed to yield the
        % covariance matrix)
        Phi_L;
        
        % Transfromed covaraince parameters
        Sigma_chol;
        
        % Log(det(Sigma))
        Log_det_Sigma;
        
        % Dimension of the distribution
        N_dim;
        
        % Total number of parameters:
        N_var_param;
    end
    
    methods
        % Constructor:
        function this = varDistGaussian(phi_mu_init, phi_L_init)
            this.N_dim = length(phi_mu_init);
            this.N_var_param = length(phi_mu_init) + length(phi_L_init);
            
            this.Phi_mu = phi_mu_init;
            this.Phi_L = phi_L_init;
            this = this.getSigmaCholFromPhi_L();
            this.Log_det_Sigma = 2 * sum(log(diag(this.Sigma_chol)));
        end
        
%         function this = setPhi(this, phi_mu, phi_L)
%             this.N_dim = length(phi_mu_init);
%             this.N_var_param = length(phi_mu_init) + length(phi_L_init);
%             this.Phi_mu = phi_mu;
%             this.Phi_L = phi_L;
%             this = this.getSigmaCholFromPhi_L();
%             this.Log_det_Sigma = 2 * sum(log(diag(this.Sigma_chol)));
%         end
        
        function this = updatePhi(this, delta_phi)
            this.Phi_mu = this.Phi_mu + delta_phi(1:this.N_dim);
            this.Phi_L = this.Phi_L + delta_phi(this.N_dim+1:end);
            this = this.getSigmaCholFromPhi_L();
            this.Log_det_Sigma = 2 * sum(log(diag(this.Sigma_chol)));
        end
        
        function X = reparamFun(this, epsilon)
            % make sure epsilon is a row vector
            if(numel(epsilon) == length(epsilon))
                epsilon = reshape(epsilon,1,numel(epsilon));
            else
               error('Incorrect format of epsilon') 
            end
            
            X = this.Phi_mu + epsilon * this.Sigma_chol';
        end
        
        function DXDphi = gradReparamFun_dphi(this, epsilon)
            DXDphi_mu = eye(this.N_dim);
            DXDPhi_L = this.getDXDPhi_L(epsilon);
            DXDphi = [DXDphi_mu, DXDPhi_L];
        end
        
        function entropy = getEntropy(this)
            entropy = 0.5 * (this.N_dim * log(2*pi*exp(1)) + this.Log_det_Sigma);
        end
        
        function DEntropyDphi = getGradEntropy_dphi(this)
            % Distinguish three different formats of Phi_L
            if(numel(this.Phi_L) == 1)
                DEntropyDphi_L = this.N_dim;
                
            elseif(numel(this.Phi_L) == this.N_dim)
                DEntropyDphi_L = ones(1,this.N_dim);
                
            elseif(numel(this.Phi_L) == this.N_dim * (this.N_dim + 1) /2)
                DEntropyDphi_L = ones(1,this.N_dim);
                % The offdiagonal terms of the covariance matrix do not
                % influence the entropy of the Gaussian.
                DEntropyDphi_L = [DEntropyDphi_L, zeros(1,length(this.Phi_L) - this.N_dim)];
                
            else
                error('Phi_L has inappropriate format')
            end
            
            DEntropyDphi_mu = zeros(1,this.N_dim);
            DEntropyDphi = [DEntropyDphi_mu, DEntropyDphi_L];
        end
        
        function this = getSigmaCholFromPhi_L(this)
            % Sigma is parametrized through Phi_L as follows:
            % Sigma = L' * L;
            %       exp(Phi_L_1)    0               0               ...
            %       Phi_L_6         exp(Phi_L_2)    0               ...
            %   L = Phi_L_7         Phi_L_8         exp(Phi_L_3)    ...
            %       Phi_L_9         Phi_L_10        Phi_L_11        ...
            %       Phi_L_12        Phi_L_13        Phi_L_14        ...
            
            if(numel(this.Sigma_chol) == 0)
                this.Sigma_chol = zeros(this.N_dim);
            end
            
            % Distinguish three different formats of Phi_L
            if(numel(this.Phi_L) == 1)
                this.Sigma_chol = eye(this.N_dim) * exp(this.Phi_L);
                
            elseif(numel(this.Phi_L) == this.N_dim)
                this.Sigma_chol = diag( exp(this.Phi_L) );
                
            elseif(numel(this.Phi_L) == this.N_dim * (this.N_dim + 1) /2)
                Sigma_chol_diag = diag( exp(this.Phi_L(1:this.N_dim)) );
                Sigma_chol_offdiag = triu(ones(this.N_dim),1);
                Sigma_chol_offdiag(Sigma_chol_offdiag == 1) = this.Phi_L(this.N_dim+1:end);
                Sigma_chol_offdiag = Sigma_chol_offdiag';
                this.Sigma_chol = Sigma_chol_diag + Sigma_chol_offdiag;
                
            else
                error('Phi_L has inappropriate format')
            end
            
        end
        
        function DXDPhi_L = getDXDPhi_L(this, epsilon)
            if(numel(this.Phi_L) == 1)
                DXDPhi_L = diag( exp(this.Phi_L) .* epsilon );
                
            elseif(numel(this.Phi_L) == this.N_dim)
                DXDPhi_L = diag( exp(this.Phi_L) .* epsilon );
                
            elseif(numel(this.Phi_L) == this.N_dim * (this.N_dim + 1) /2)
                DXDPhi_L = diag( exp(this.Phi_L(1:this.N_dim)) .* epsilon );
                
                help_mat = triu(ones(this.N_dim),1);
                help_mat(help_mat == 1) = 1:(this.N_dim-1) * (this.N_dim) /2;
                help_mat = help_mat';
                help_mat(help_mat>0) = help_mat(help_mat>0) + this.N_dim;
                
                for row = 2:this.N_dim
                    indxs_Phi = help_mat(row,help_mat(row,:)>0);
                    indxs_eps = 1:row-1;
                    DXDPhi_L(row, indxs_Phi) = epsilon(indxs_eps);
                end
                
            else
                error('Phi_L has inappropriate format')
            end
        end
        
        function X_samples = generateSamples(this, N_samples)
            % Generate samples from the standard normal distribution:
            standard_samples = randn(N_samples, this.N_dim);
            
            % transfrom the standard samples
            X_samples = repmat(this.Phi_mu, N_samples, 1)...
                + standard_samples * this.Sigma_chol';
        end
        
    end
    
end