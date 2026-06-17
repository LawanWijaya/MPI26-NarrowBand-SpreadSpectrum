function return_data = test_basin(bestfit)

    gap = bestfit-[0,bestfit(1:end-1)];
    gap = gap(2:end);
    half_gap = gap/2;

    basin = zeros(1,length(bestfit));

    basin(1) = bestfit(1) + half_gap(1);
    for i = 2:length(bestfit)-1
        basin(i) = half_gap(i-1) + half_gap(i);
    end
    basin(end) = half_gap(end) + max(bestfit);

    return_data = min(basin);

end