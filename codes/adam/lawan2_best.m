function [selectedIdx, selectedFreq, bestGap, gapVec] = lawan2_best(faithfulIdx, w)
% Select K faithful Walsh rows so that:
%   1. the rows are spread apart in Walsh frequency/sign-change count,
%   2. each consecutive frequency gap is odd, and
%   3. the minimum consecutive gap is as large as possible.

K = 2^w;

if length(faithfulIdx) < K
    error('Not enough faithful rows. Need %d rows, but only %d are faithful.', ...
        K, length(faithfulIdx));
end

% faithfulFreq = rowFreq(faithfulIdx);
faithfulFreq = faithfulIdx-1;


[faithfulFreqSorted, order] = sort(faithfulFreq);
faithfulIdxSorted = faithfulIdx(order);

low = 0;
high = max(faithfulFreqSorted) - min(faithfulFreqSorted);

bestGap = 0;
bestSelectedIdx = [];

while low <= high
    gap = floor((low + high)/2);

    % Now the greedy selection also enforces that every consecutive
    % frequency difference is odd.
    chosenIdx = greedyChooseWithOddGap(faithfulIdxSorted, faithfulFreqSorted, K, gap);

    if length(chosenIdx) >= K
        bestSelectedIdx = chosenIdx(1:K);
        low = gap + 1;
    else
        high = gap - 1;
    end
end

if isempty(bestSelectedIdx)
    error('Could not find %d faithful rows with odd consecutive frequency gaps.', K);
end

selectedIdx = bestSelectedIdx(:);
% selectedFreq = rowFreq(selectedIdx);
selectedFreq = selectedIdx-1;

% Sort the final answer by frequency so that the displayed gaps are consecutive.
[selectedFreq, order2] = sort(selectedFreq);
selectedIdx = selectedIdx(order2);

gapVec = diff(selectedFreq); % Need the gap between selected frequencies to see if they are odd

% Safety check: every consecutive frequency gap should be odd.
if any(mod(gapVec,2) == 0)
    error('Odd-gap constraint failed: at least one consecutive frequency gap is even.');
end

% Report the actual minimum consecutive odd gap for the selected set.
bestGap = min(gapVec);

end


function chosenIdx = greedyChooseWithOddGap(idxSorted, freqSorted, K, gap)
% Greedily choose rows with:
%   frequency difference >= gap
%   frequency difference odd
%
% We try every possible starting row. This avoids forcing the selection to
% start from the very first faithful row, which can be too restrictive once
% the odd-gap constraint is added.

bestChosenIdx = [];

for start = 1:length(idxSorted)

    chosenIdx = idxSorted(start);
    lastFreq = freqSorted(start);

    for j = start+1:length(idxSorted)

        actualGap = freqSorted(j) - lastFreq;

        if actualGap >= gap && mod(actualGap,2) == 1
            chosenIdx(end+1) = idxSorted(j); %#ok<AGROW>
            lastFreq = freqSorted(j);

            if length(chosenIdx) == K
                return;
            end
        end
    end

    if length(chosenIdx) > length(bestChosenIdx)
        bestChosenIdx = chosenIdx;
    end
end

chosenIdx = bestChosenIdx;

end
