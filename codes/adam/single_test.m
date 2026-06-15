function return_data = single_test(sz,cutoff)

    inp_seq = 2*randi([0,1],sz,1) - 1;
    out_seq = adam_faithful(inp_seq,cutoff);
    diff_counts = sum(inp_seq ~= out_seq, 1);

    return_data = diff_counts;

end