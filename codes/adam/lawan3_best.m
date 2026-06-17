function selectedIdx = lawan3_best(faithfulIdx, w)

K = 2^w;

if length(faithfulIdx) < K
    error('Not enough faithful rows. Need %d rows, but only %d are faithful.', ...
        K, length(faithfulIdx));
end

faithfulFreq = faithfulIdx-1;

[faithfulFreqSorted, order] = sort(faithfulFreq);
faithfulIdxSorted = faithfulIdx(order);

low = 0;
high = max(faithfulFreqSorted) - min(faithfulFreqSorted);
bestGap = 0;
bestSelectedIdx = [];

while low <= high
    gap = floor((low + high)/2);

    chosenIdx = greedyChooseWithGap(faithfulIdxSorted, faithfulFreqSorted, K, gap);

    if length(chosenIdx) >= K
        bestGap = gap;
        bestSelectedIdx = chosenIdx(1:K);
        low = gap + 1;
    else
        high = gap - 1;
    end
end

selectedIdx = bestSelectedIdx(:);
% selectedFreq = rowFreq(selectedIdx);
selectedFreq = selectedIdx-1;


[selectedFreq, order2] = sort(selectedFreq);
selectedIdx = selectedIdx(order2);
end


function chosenIdx = greedyChooseWithGap(idxSorted, freqSorted, K, gap)

chosenIdx = idxSorted(1);
lastFreq = freqSorted(1);

for j = 2:length(idxSorted)
    if freqSorted(j) - lastFreq >= gap
        chosenIdx(end+1) = idxSorted(j); %#ok<AGROW>
        lastFreq = freqSorted(j);

        if length(chosenIdx) == K
            return;
        end
    end
end
end