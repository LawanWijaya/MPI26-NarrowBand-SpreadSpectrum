function return_data = simulate_fading(inpseq,options)
% SIMULATE_FADING functions call take two forms: either submit only an
% input sequence, and default probabilities will be applied, or manually
% specify the probabilities.
arguments (Input)
    inpseq                  % the input sequence
end
arguments (Input)
    options.stay_good = 0.9  % probability of staying faithful
    options.stay_bad  = 0.5  % probability of propagating error
end

    state = 0;
    outseq = inpseq;

    for i = 1:length(inpseq)
        if state % preceding i has error
            if rand > options.stay_bad   % error ends
                state = 0;
                disp(i)
            else                         % error persists
                outseq(i) = -1*inpseq(i);
            end
        else     % preceding i is correct
            if rand > options.stay_good  % begin error
                disp(i)
                state = 1;
                outseq(i) = -1*inpseq(i);
            end
        end
    end

    return_data = outseq;

end