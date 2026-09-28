function bestAlignCode = score_all_alignments()
clear; clc;

fprintf("Script started\n");

if ~exist('outputs', 'dir')
    mkdir('outputs');
end

binSize = 0.01;
window = [-0.5, 1];
edges = window(1):binSize:window(2);
alignCodes = [0, 1, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 44, 45, 46, 48, 49, ...
              50, 51, 52, 53, 54, 55, 56, 57, 64, 65, 66, 255];

files = dir('data/**/*.mat');
files = files(~contains({files.name}, 'plxAnaAD'));

results = [];

for a = 1:length(alignCodes)
    alignCode = alignCodes(a);
    allCounts = [];

    fprintf("Testing alignCode = %d...\n", alignCode);

    for k = 1:length(files)
        try
            s = load(fullfile(files(k).folder, files(k).name));
            if isfield(s, 'data') && isfield(s.data, 'chan1a') && isfield(s, 'strobed')
                spikes = s.data.chan1a;
                strobed = s.strobed;

                alignTimes = strobed(strobed(:,2) == alignCode, 1);
                if isempty(alignTimes), continue; end

                [counts, ~] = computePSTH(spikes, alignTimes, edges);
                allCounts = [allCounts; counts];
            end
        catch
            continue;
        end
    end

    if isempty(allCounts)
        fprintf("No usable data for alignCode = %d\n", alignCode);
        continue;
    end

    fr = mean(allCounts, 1) / binSize;
    binCenters = edges(1:end-1) + diff(edges(1:2))/2;

    peakFR = max(fr);
    peakIdx = find(fr == peakFR, 1);
    peakTime = binCenters(peakIdx);
    baseline = mean(fr(binCenters < -0.2));
    stdBL = std(fr(binCenters < -0.2));
    if stdBL == 0, stdBL = 1e-6; 
    end
    skewnessScore = mean(fr(binCenters > 0)) - mean(fr(binCenters < 0));

    score = (peakFR - baseline) / stdBL + (abs(peakTime) < 0.1)*2 + skewnessScore;

    results = [results; alignCode, score, peakFR, peakTime];
end

% Heuristic ranking only; selection on these data is exploratory.
if isempty(results)
    error('No usable events found. Check authorized input data and event codes.');
end

% Sort and save
[~, idx] = sort(results(:,2), 'descend');
results_sorted = results(idx, :);

T = array2table(results_sorted, 'VariableNames', {'AlignCode', 'Score', 'PeakFR', 'PeakTime'});
writetable(T, 'outputs/alignment_scores.csv');

bestAlignCode = results_sorted(1,1);
fprintf("Scoring complete. Best alignCode = %d\n", bestAlignCode);
end
