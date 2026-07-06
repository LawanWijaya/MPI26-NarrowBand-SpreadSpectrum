function return_data = test_basin(bestfit)
% TEST_BASIN computes the size of the basin for each element in a vector.
% A 'basin' is defined as the set of positive integers closer to a given
% than any other.
%
%last updated 06/26/2026 by Adam Petrucci
arguments
    bestfit       % the vector to compute basins for
end

    % Compute differences between subsequent elements in the vector
    gap = bestfit-[0,bestfit(1:end-1)];
    gap = gap(2:end);
    half_gap = gap/2;

    % Initalize storage objects
    basin = zeros(1,length(bestfit));

    % Compute basins for each element of vector
    basin(1) = bestfit(1) + half_gap(1); % first element 'attracts'
                                         % everything preceding it
    for i = 2:length(bestfit)-1
        basin(i) = half_gap(i-1) + half_gap(i); % interior elements
                                                % attract the nearest half
                                                % of their adjacent gaps
    end
    basin(end) = half_gap(end) + max(bestfit); % last element 'attracts'
                                               % everything following it
    
    % Note: by design of the walsh ordering, the last element always has
    % the largest basin, so is not of interest when computing the min

    % Return minimum basin size
    return_data = min(basin);

end