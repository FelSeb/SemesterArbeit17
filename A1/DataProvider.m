classdef DataProvider
    %DATAPROVIDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Data directory (for now we store the training data directly)
        Data_directory;
        I2B;
        B2I;
        current_batch;
        
    end
    
    methods
        % Constructor
        function this = DataProvider(data_directory)
            this.Data_directory = data_directory;
        end
        
        % This funciton returns training sample i. It loads training data
        % batchwise into the workspace in order to economize storage space.
        function [x_i, y_i] = provideDataPoint(this,i)
            
            % Load index to batch map and its inverse
            if(isempty(this.I2B))
                load([this.Data_directory,'/I2B']);
                load([this.Data_directory,'/B2I']);
                this.I2B = I2B;
                this.B2I = B2I;
            end
            
            % load new batch
            if(isempty(this.current_batch) || i == 1 || this.I2B(1,i) > this.I2B(1,i-1))
                load([this.Data_directory,'/data_point_batch',num2str(this.I2B(1,i))])
                this.current_batch = data_point_batch;
            end
            % extract one data point from batch
            x_i = this.current_batch{this.I2B(2,i)}.x;
            y_i = this.current_batch{this.I2B(2,i)}.y;
        end
    end
    
end

