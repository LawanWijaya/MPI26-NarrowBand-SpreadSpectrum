function return_data = david2_best(feasible,w)
% DAVID2_BEST generates a codebook according to David's first algorithm
%
%last updated 06/28/2026 by Adam Petrucci

    greedy = david_best(feasible,w);

    return_data = greedy(:,2);

end