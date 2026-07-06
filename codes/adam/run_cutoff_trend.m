function [return_data_coarse, return_data_fine] = ...
                                               run_cutoff_trend(sz,cutoffs)
% RUN_CUTOFF_TREND computes the faithful rows of a Hadamard matrix (with
% Walsh ordering) given the cutoff associated to a low-pass filter.
%
%last updated 06/26/2026 by Adam Petrucci
arguments
    sz         % length of encoding chip sequence
    cutoffs    % truncation values
end

    % Construct Walsh matrix (Hadamard rows in Walsh order)
    N = 2^sz;
    Hn = HadtoW(sz);

    % Initialize storage objects
    fcs = zeros(1, length(cutoffs));        % number of faithful rows
    fine_result = zeros(N,length(cutoffs)); % faithful rows

    % Determine faithful set for all given cutoff values
    for c_e = 1:length(cutoffs)

        % Select cutoff and reset counter
        cf = cutoffs(c_e);
        faithful_count = 0;

        % Test each row in the matrix
        for h = 1:N

            % Select row (inp_seq = input sequence)
            % and apply low-pass filter (out_seq = output sequence)
            inp_seq = Hn(h,:).';
            out_seq = adam_faithful(inp_seq,cf);

            % Cross-correlate ouput sequence with rows of Hadamard by
            % computing boolean distance
            diff_counts = sum(Hn ~= out_seq, 1);

            % Find rows of minimum distance from output sequence
            % (min_val = minimum value, idxs = indices)
            min_val = min(diff_counts);
            idxs = find(diff_counts == min_val);

            % Check if 'optimal' row (by cross-correlation) is input
            % sequence (close_match = true/false for affirmative/negative)
            close_match = 0;
            if isscalar(idxs)    % first check that 'optimal' row is unique
                if inp_seq == Hn(idxs,:).' % then compare to input sequence
                    close_match = 1;
                end
            end

            % If input sequence was preserved, update storage objects
            if close_match
                faithful_count = faithful_count + 1; % increment count of
                                                     % faithful rows
                fine_result(h,c_e) = 1; % record faithful row
            end

        end

        % Record total number of faithful rows
        fcs(c_e) = faithful_count;

    end

    % Return all data
    return_data_coarse = fcs;
    return_data_fine = fine_result;

end