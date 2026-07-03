function return_data = david1_best(feasible,w)
% DAVID1_BEST generates a codebook according to David's first algorithm
%
%last updated 06/28/2026 by Adam Petrucci

    greedy = david_best(feasible,w);

    return_data = greedy(:,1);

end