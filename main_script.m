
% Event-aligned exploratory analysis; condition labels require validation.
% Merge raster and PSTH across all usable .mat files, with vision comparison

clear; clc;

% Create output directory if it doesn't exist
if ~exist('outputs', 'dir')
    mkdir('outputs')
end




% Parameters
binSize = 0.01;
window = [-0.5, 1];
edges = window(1):binSize:window(2);
% Configured event code; verify its meaning against acquisition metadata.
% score_all_alignments() provides an optional exploratory ranking.
% No saved ranking is included to establish why code 0 was selected.
alignCode = 0;

files = dir('data/**/*.mat');
files = files(~contains({files.name}, 'plxAnaAD'));
if isempty(files)
    error('No input files found. Supply authorized MAT files under data/. See README.');
end

allSpikeCounts = [];
visionSpikeCounts = [];
noVisionSpikeCounts = [];
allSpikesPerTrial = {};
visionSpikesPerTrial = {};
noVisionSpikesPerTrial = {};

trialCounter = 0;
visionCounter = 0;
noVisionCounter = 0;

for k = 1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    fprintf("Checking %s\n", filePath);
    try
        s = load(filePath);
        if isfield(s, 'data') && isfield(s.data, 'chan1a') && isfield(s, 'strobed')
            spikes = s.data.chan1a;
            strobed = s.strobed;

            alignTimes = strobed(strobed(:,2) == alignCode, 1);
            if isempty(alignTimes)
                fprintf("No alignCode=%d in %s\n", alignCode, files(k).name);
                continue;
            end

            % ASSUMPTION: odd/even events represent Vision/No Vision.
            % Verify against experimental metadata before interpreting conditions.
            visionAlignTimes = alignTimes(1:2:end);
            noVisionAlignTimes = alignTimes(2:2:end);

            % Collect spike data
            for i = 1:length(alignTimes)
                aligned = spikes - alignTimes(i);
                aligned = aligned(aligned >= window(1) & aligned <= window(2));
                allSpikesPerTrial{end+1} = aligned;
            end

            for i = 1:length(visionAlignTimes)
                aligned = spikes - visionAlignTimes(i);
                aligned = aligned(aligned >= window(1) & aligned <= window(2));
                visionSpikesPerTrial{end+1} = aligned;
            end

            for i = 1:length(noVisionAlignTimes)
                aligned = spikes - noVisionAlignTimes(i);
                aligned = aligned(aligned >= window(1) & aligned <= window(2));
                noVisionSpikesPerTrial{end+1} = aligned;
            end

            [c_all, ~] = computePSTH(spikes, alignTimes, edges);
            allSpikeCounts = [allSpikeCounts; c_all];

            [c_v, ~] = computePSTH(spikes, visionAlignTimes, edges);
            visionSpikeCounts = [visionSpikeCounts; c_v];

            [c_nv, ~] = computePSTH(spikes, noVisionAlignTimes, edges);
            noVisionSpikeCounts = [noVisionSpikeCounts; c_nv];

        else
            fprintf("Skipped: missing fields\n");
        end
    catch ME
        fprintf("Error in %s: %s\n", files(k).name, ME.message);
    end
end

% Diagnostic: how many non-empty trials
nonemptyTrials = sum(cellfun(@(x) ~isempty(x), allSpikesPerTrial));
fprintf('Total non-empty raster trials: %d\n', nonemptyTrials);

% Plot 1: merged raster
if nonemptyTrials > 0
    fig1 = figure; hold on;
    for t = 1:length(allSpikesPerTrial)
        ts = allSpikesPerTrial{t};
        if ~isempty(ts)
            plot(ts, t*ones(size(ts)), 'k.');
        end
    end
    title('Merged Raster');
    xlabel('Time (s)');
    ylabel('Trial');
    xlim(window);
    saveas(fig1, 'outputs/merged_raster.png');
    close(fig1);
end

% Plot 2: overall PSTH
if ~isempty(allSpikeCounts)
    totalPSTH = mean(allSpikeCounts,1)/binSize;
    binCenters = edges(1:end-1) + diff(edges(1:2))/2;
    fig2 = figure;
    bar(binCenters, totalPSTH, 'FaceColor', 'k');
    title('Merged PSTH');
    xlabel('Time (s)');
    ylabel('Firing Rate (Hz)');
    saveas(fig2, 'outputs/merged_psth.png');
    close(fig2);
end

% Plot 3: Vision vs No Vision PSTH
if ~isempty(visionSpikeCounts) && ~isempty(noVisionSpikeCounts)
    visionPSTH = mean(visionSpikeCounts,1)/binSize;
    noVisionPSTH = mean(noVisionSpikeCounts,1)/binSize;
    binCenters = edges(1:end-1) + diff(edges(1:2))/2;

    % Difference threshold before smoothing
    diffPSTH = abs(visionPSTH - noVisionPSTH);
    threshold = 5;
    differenceBins = diffPSTH > threshold;

    % Smoothing
    visionPSTH = smoothdata(visionPSTH, 'gaussian', 5);
    noVisionPSTH = smoothdata(noVisionPSTH, 'gaussian', 5);

    fig3 = figure; hold on;
    plot(binCenters, visionPSTH, 'b', 'LineWidth', 2);
    plot(binCenters, noVisionPSTH, 'r', 'LineWidth', 2);
    title('Vision vs No Vision PSTH');
    xlabel('Time (s)');
    ylabel('Firing Rate (Hz)');


    % Descriptive markers: not a statistical test or significance claim.
    markerY = max([visionPSTH, noVisionPSTH]) * 1.1;
    plot(binCenters(differenceBins), markerY * ones(1, sum(differenceBins)), 'k.');

    legend('Vision (assumed odd events)', 'No Vision (assumed even events)', ...
        'Absolute unsmoothed difference > 5 Hz (descriptive)');

    % Save descriptive summaries; no p-values or inferential test.
    stats = struct();
    stats.binCenters = binCenters;
    stats.visionPSTH = visionPSTH;
    stats.noVisionPSTH = noVisionPSTH;
    stats.diffPSTH = diffPSTH;
    stats.differenceBins = differenceBins;
    stats.threshold = threshold;

    [stats.visionPeakFR, idxV] = max(visionPSTH);
    [stats.noVisionPeakFR, idxNV] = max(noVisionPSTH);
    stats.visionPeakTime = binCenters(idxV);
    stats.noVisionPeakTime = binCenters(idxNV);

    save('outputs/vision_comparison_summary.mat', 'stats');
    fprintf('Saved descriptive summary to outputs/vision_comparison_summary.mat\n');

    saveas(fig3, 'outputs/vision_comparison_descriptive.png');
    close(fig3);
end
