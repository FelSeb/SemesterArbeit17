clear;
clc;

nsamples = 10;
N_batch = 11;

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

