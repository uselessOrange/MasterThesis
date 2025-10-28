function results = analyzeCompressorSegments_legacy(t, x, compressor, usableDataIndex, varargin)
% analyzeCompressorSegments
%   Performs peak/valley analysis, compressor ON/OFF mask prediction,
%   and comparison against measured compressor data for multiple segments.
%
%   results = analyzeCompressorSegments(t, x, compressor, usableDataIndex)
%   x-Temperature in the fridge in time

%   t-Time vactor

%   compressor-Recorded compressor states to be compared with predicted
%   states

%   usableDataIndex-Index of values to be used in prediction after anomaly
%   exclusion


% Optional Name-Value arguments:
%   'MinPeakProminence' — prominence threshold for findpeaks (default 0.8)
%   'MinPeakDistance'   — minimum distance between peaks (default 0.17)
%   'PlotResults'       — true/false, show figures (default true)
%
% Output: results — struct array with fields:
%   .pks, .locs_peaks, .valleys, .locs_valleys
%   .mask_pred, .mask_meas
%   .stats (struct from compareOnTime)
%   .t_chunk, .x_chunk, .compressor_chunk

%% --- Parse inputs ---
p = inputParser;
addParameter(p, 'MinPeakProminence', 0.8, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'MinPeakDistance',   0.17, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'PlotResults',       true, @(x)islogical(x)||ismember(x,[0,1]));
parse(p, varargin{:});
params = p.Results;

%% --- Initialization ---
nSets = numel(usableDataIndex);
results = struct('pks', [], 'locs_peaks', [], 'valleys', [], 'locs_valleys', [], ...
                 'mask_pred', [], 'mask_meas', [], 'stats', [], ...
                 't_chunk', [], 'x_chunk', [], 'compressor_chunk', []);

%% --- Loop through each segment ---
for i = 1:nSets
    idx = usableDataIndex{i};
    N = numel(idx);

    % Extract current chunk
    t_chunk = t(idx);
    x_chunk = x(idx);
    compressor_chunk = compressor(idx);

    % --- 1) Find peaks and valleys ---
    [pks, locs_peaks, valleys, locs_valleys] = ...
        findPeaksAndValleys(x_chunk, ...
            'MinPeakProminence', params.MinPeakProminence, ...
            'MinPeakDistance',   params.MinPeakDistance);

    % --- 2) Build predicted mask ---
    mask = peakValleyMask(locs_peaks, locs_valleys, N);

    % --- 3) Align from first peak to end ---
    startPK = locs_peaks(1);
    mask_pred = mask(startPK:end);
    mask_meas = compressor_chunk(startPK:end);
    mask_meas = mask_meas(1:numel(mask_pred));

    % --- 4) Compare predicted vs measured ---
    stats = compareOnTime(mask_pred, mask_meas);

    % --- 5) Save results ---
    results(i).pks            = pks;
    results(i).locs_peaks     = locs_peaks;
    results(i).valleys        = valleys;
    results(i).locs_valleys   = locs_valleys;
    results(i).mask_pred      = mask_pred;
    results(i).mask_meas      = mask_meas;
    results(i).stats          = stats;
    results(i).t_chunk        = t_chunk;
    results(i).x_chunk        = x_chunk;
    results(i).compressor_chunk = compressor_chunk;

    % --- 6) Plot results if requested ---
    if params.PlotResults
        % (a) Peak and valley plot
        figure('Name', sprintf('Peaks/Valleys Segment %d', i))
        plot(t_chunk/3600, x_chunk, 'b-'); hold on
        plot(t_chunk(locs_peaks)/3600, pks, 'ro', 'MarkerFaceColor', 'r');
        plot(t_chunk(locs_valleys)/3600, valleys, 'go', 'MarkerFaceColor', 'g');
        xlabel('Time (h)'), ylabel('Amplitude')
        title(sprintf('Peak Detection — Segment %d', i))
        grid on, hold off

        % (b) Compressor ON/OFF comparison
        figure('Name', sprintf('Compressor Compare Segment %d', i))
        plot(mask_pred, '-r', 'LineWidth', 2); hold on
        plot(mask_meas, '-b'); hold off
        xlabel('Sample index'), ylabel('Compressor ON/OFF')
        legend('Predicted','Measured')
        title(sprintf('Compressor ON/OFF Comparison — Segment %d', i))
        grid on
    end

    % --- 7) Print stats summary ---
    fprintf('--- Segment %d ---\n', i)
    fprintf('Measured ON-time   = %.1f h\n', stats.timeOnMeas/3600)
    fprintf('Predicted ON-time  = %.1f h\n', stats.timeOnPred/3600)
    fprintf('Difference         = %+0.1f s (%.1f%%)\n\n', ...
            stats.diffTime, stats.pctError)
end

end