function [counts, binEdges] = computePSTH(spikes, alignTimes, edges)
    % Optimized PSTH computation
    numTrials = length(alignTimes);
    numBins = length(edges) - 1;
    counts = zeros(numTrials, numBins);
    
    % Use vectorized computation when possible
    for i = 1:numTrials
        % Compute aligned spike times once
        alignedSpikes = spikes - alignTimes(i);
        % Use builtin histcounts for speed
        counts(i, :) = histcounts(alignedSpikes, edges);
    end
    
    binEdges = edges;
end