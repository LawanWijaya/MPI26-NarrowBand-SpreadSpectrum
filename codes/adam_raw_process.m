addpath(genpath('C:\Users\adamn\OneDrive\Documents\GitHub\MPI26-NarrowBand-SpreadSpectrum'))


%%

inp_seq = 2*randi([0,1],8,1) - 1;
out_seq = adam_faithful(inp_seq,1);
num_diff = sum(inp_seq ~= out_seq)