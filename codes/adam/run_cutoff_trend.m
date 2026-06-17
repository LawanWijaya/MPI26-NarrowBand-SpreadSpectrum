function [return_data_coarse, return_data_fine] = run_cutoff_trend(sz,cutoffs)

    N = 2^sz;
    %Hn = hadamard(N);
    Hn = HadtoW(sz);
    %Hn = HadtoD(sz);

    fcs = zeros(1, length(cutoffs));

    fine_result = zeros(sz,length(cutoffs));

    for c_e = 1:length(cutoffs)

        cf = cutoffs(c_e);

        faithful_count = 0;

        for h = 1:N

            inp_seq = Hn(h,:).';
            out_seq = adam_faithful(inp_seq,cf);

            diff_counts = sum(Hn ~= out_seq, 1);

            min_val = min(diff_counts);
            idxs = find(diff_counts == min_val);

            close_match = 0;
            if isscalar(idxs)
                if inp_seq == Hn(idxs,:).'
                    close_match = 1;
                end
            end

            if close_match
                faithful_count = faithful_count + 1;
                fine_result(h,c_e) = 1;
            end

        end

        fcs(c_e) = faithful_count;

    end

    return_data_coarse = fcs;
    return_data_fine = fine_result;

end