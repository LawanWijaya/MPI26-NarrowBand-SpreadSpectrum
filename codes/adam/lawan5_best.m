function selectedIdx = lawan5_best(faithfulIdx, w)
% LAWAN5_BEST generates a codebook from given a given feasible set.
% Everything below this comment block was emailed to Adam from Lawan on
% 6/17/26 and has not been modified.
%
%last updated 06/28/2026 by Adam Petrucci

K = 2^w;

if length(faithfulIdx) < K
    error('Not enough faithful rows. Need %d rows, but only %d are faithful.', ...
        K, length(faithfulIdx));
end

% Frequency / sequency of each faithful Walsh row
faithfulFreq = faithfulIdx - 1;

% Sort faithful rows by increasing frequency
[faithfulFreqSorted, order] = sort(faithfulFreq);
faithfulIdxSorted = faithfulIdx(order);

% Binary search for the largest possible minimum frequency gap
low = 0;
high = max(faithfulFreqSorted) - min(faithfulFreqSorted);

bestGap = 0;
bestSelectedIdx = [];

while low <= high

    gap = floor((low + high)/2);

    % Try to select K rows starting from the maximum frequency
    % and moving backward
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

% Sort selected rows in increasing frequency order before returning
% This makes the output easier to read and consistent with the codebook format
selectedFreq = selectedIdx - 1;
[selectedFreq, order2] = sort(selectedFreq);
selectedIdx = selectedIdx(order2);

end


function chosenIdx = greedyChooseWithGap(idxSorted, freqSorted, K, gap)

% Start from the highest-frequency faithful row
chosenIdx = idxSorted(end);
lastFreq = freqSorted(end);

% Move backward toward lower frequencies
for j = length(idxSorted)-1:-1:1

    if lastFreq - freqSorted(j) >= gap

        chosenIdx(end+1) = idxSorted(j); %#ok<AGROW>
        lastFreq = freqSorted(j);

        if length(chosenIdx) == K
            return;
        end

    end

end

end