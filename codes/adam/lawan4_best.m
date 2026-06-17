function selectedIdx = lawan4_best(faithfulIdx, w)

    K = 2^w;

    if length(faithfulIdx) < K
        error('Not enough faithful rows. Need %d rows, but only %d are faithful.', ...
            K, length(faithfulIdx));
    end

    % Frequency is assumed to be row index minus 1
    faithfulFreq = faithfulIdx - 1;

    % Sort by frequency
    [faithfulFreqSorted, order] = sort(faithfulFreq);
    faithfulIdxSorted = faithfulIdx(order);

    low = 0;
    high = max(faithfulFreqSorted) - min(faithfulFreqSorted);

    bestGap = 0;
    bestSelectedIdx = [];

    while low <= high

        gap = floor((low + high)/2);

        chosenIdx = dynamicChooseWithGap( ...
            faithfulIdxSorted, faithfulFreqSorted, K, gap);

        if length(chosenIdx) == K
            bestGap = gap;
            bestSelectedIdx = chosenIdx;
            low = gap + 1;
        else
            high = gap - 1;
        end

    end

    if isempty(bestSelectedIdx)
        error('Could not find a valid set of %d rows.', K);
    end

    selectedIdx = bestSelectedIdx(:);

    % Sort final answer by frequency
    selectedFreq = selectedIdx - 1;
    [~, order2] = sort(selectedFreq);
    selectedIdx = selectedIdx(order2);

    %fprintf('Best minimum frequency gap = %d\n', bestGap);

end


function chosenIdx = dynamicChooseWithGap(idxSorted, freqSorted, K, gap)

    m = length(idxSorted);

    % dp(j,t) = true if we can choose t rows ending at row j
    dp = false(m, K);

    % parent(j,t) stores the previous row index before j
    parent = zeros(m, K);

    % Any row can be the first selected row
    dp(:,1) = true;

    for t = 2:K

        for j = 1:m

            % Try all previous rows i before j
            for i = 1:j-1

                if dp(i,t-1) && freqSorted(j) - freqSorted(i) >= gap

                    dp(j,t) = true;
                    parent(j,t) = i;

                    % Stop once we find one valid predecessor
                    break;

                end

            end

        end

    end

    % Check whether any row can end a valid K-row sequence
    endCandidates = find(dp(:,K));

    if isempty(endCandidates)
        chosenIdx = [];
        return;
    end

    % Pick one valid ending row.
    % This does not force the sequence to start from the lowest frequency.
    endRow = endCandidates(1);

    % Reconstruct selected rows
    selectedRows = zeros(K,1);
    selectedRows(K) = endRow;

    for t = K:-1:2
        selectedRows(t-1) = parent(selectedRows(t), t);
    end

    chosenIdx = idxSorted(selectedRows);

end