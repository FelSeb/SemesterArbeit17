classdef DataBatch
    %DATABATCH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Batch_size;
        Data_point_batch;
    end
    
    methods
        function  this = DataBatch(batch_size)
           this.Batch_size = batch_size;
           this.Data_point_batch = cell(1,batch_size);
        end
        
        function this = addDataPoint_i(this, DP,i)
           this.Data_point_batch{i} = DP;
        end
        
    end
    
end

