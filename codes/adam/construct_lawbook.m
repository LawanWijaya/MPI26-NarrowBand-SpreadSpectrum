function return_data = construct_lawbook(w,faithful,find_best)
% CONSTRUCT_LAWBOOK greedily selects a maximal set of codebooks from a set
% of faithful sequences using a given selection algorithm. 'Maximal' in
% this setting means as many distinct (disjoint) codebooks as can be
% constructed from the given faithful set with the constraint on codebook
% size
%
%last updated 06/26/2026 by Adam Petrucci
arguments
    w           % length of word (ie bits/word)
    faithful    % set of faithful rows
    find_best   % algorithm for constructing codebook; algorithm must take
                % two inputs: the set of faithful rows and the word size,
                % and must return one output: a vector of walsh frequencies
end

    % Maximum number of codebooks of size 2^w with given faithful set
    % (cds = codebooks)
    cds = floor(length(faithful)/(2^w));

    % Initialize storage objects
    lawbook = zeros(cds,2^w);

    % Greedily collect codebooks by iteratively finding best codebook then
    % removing those used sequences from the admissible set
    for cd = 1:cds

        best = find_best(faithful,w);
        lawbook(cd,:) = best;
        faithful = faithful(~ismember(faithful,best));

    end

    % Return data
    return_data = lawbook;

end