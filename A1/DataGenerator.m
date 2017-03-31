classdef DataGenerator
    %DATAGENERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        Target_dir
        
        % x and y nodal coordinates of the high resolution finite element
        % mesh. Format: Matricies with as many rows as nodes per element
        % and as many columns as elements (nodes_per_elem x N_elem x dim).
        Fine_grid;
        
        % Solver object
        Solver;
    end
    
    methods
        % Constructor
        function this = DataGenerator(fine_grid, modelSolver,dir)
            this.Fine_grid = fine_grid;
            this.Solver = modelSolver;
            % Set_up solver:
            % 1.) Set grid
            this.Solver = this.Solver.setGrid(this.Fine_grid);
            % 2.) Preevaluate shape- /test- functions (and derivatives)
            this.Solver = this.Solver.preEval();
            
            % Store
            if(~isdir(dir))
                mkdir(dir);
            end
            this.Target_dir = dir;
        end
        
        % Generate Data
        function generate(this, n_data_points, N_batch)
            % Generate
            if(n_data_points < N_batch)
               N_batch =  n_data_points;
               warning('Number of batches was reduced')
            end
            
            % Make Batch to Index and its inverse map (there are better ways)
            B2I = this.makeB2I(n_data_points, N_batch);
            I2B = this.makeI2B(B2I,n_data_points);
            save([this.Target_dir,'/B2I'],'B2I');
            save([this.Target_dir,'/I2B'],'I2B');
            
            % Iterate over batches
            for b = 1:N_batch
                if(b > 1)
                    clear data_point_batch;
                end
                
                %db = DataBatch(N_batch);
                data_point_batch = cell(1,N_batch);
                % Iterate over data points within batch
                for i = 1:length(B2I{b})
                    % Evaluate the model. A new random sample of the
                    % parameter field is drawn automatically if the
                    % parameter field is set to be a random one.
                    [output, this.Solver] = this.Solver.evaluateModel;
                    
                    % Make a DataPoint object
                    data_point_batch{i} = DataPoint(this.Solver);
                    %db = db.addDataPoint_i(data_point_batch{i}, i);
                    fprintf('Generation of datapoint no. %d complete.\n', B2I{b}(i))
                end
                
                savedir = [this.Target_dir,'/data_point_batch',num2str(b),'.mat'];
                save(savedir,'data_point_batch');
            end
        end
        
    end
    
    
    methods (Static)
        
        function B2I = makeB2I(nsamples, N_batch)
            % make the batch to index map
            from = 1;
            samp_per_batch = ceil(nsamples/N_batch);
            to = samp_per_batch;
            for b = 1:N_batch
                B2I{b} = [from:to];
                if(to+1 <= nsamples)
                    from = to+1;
                else
                    break;
                end
                
                if(to+samp_per_batch <= nsamples)
                    to = to+samp_per_batch;
                else
                    to = nsamples;
                end
            end
        end
        
        
        function I2B = makeI2B(B2I,nsamples)
            % make the index to batch map
            I2B = zeros(2,nsamples);
            for b = 1:numel(B2I)
                I2B(1,B2I{b}) = b;
                I2B(2,B2I{b}) = 1:length(B2I{b});
            end
        end
        
        
    end
end

% General To Do
% 5.) find a way to deal with data in batches in getSuffStats
% -> logQopt_i ...there should be a map i -> number_of_batch ->
% number_of_sample within batch. OR...
% Add more feature functions to pc!

%