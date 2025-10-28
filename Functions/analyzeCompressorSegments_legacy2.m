function results = analyzeCompressorSegments(t, x, compressor, usableDataIndex, varargin)
% analyzeCompressorSegments
%   Performs peak/valley analysis, compressor ON/OFF mask prediction,
%   and comparison against measured compressor data for multiple segments.
%
%   Now supports optional peak and valley index offset arguments:
%       'PeakOffset'   — integer offset to shift all peak locations
%       'ValleyOffset' — integer offset to shift all valley locations
%
% Example:
%   results = analyzeCompressorSegments(t, x, compressor, usableDataIndex, ...
%       'PeakOffset', 5, 'ValleyOffset', -3);

%% --- Parse inputs ---
p = inputParser;
addParameter(p, 'MinPeakProminence', 0.8, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'MinPeakDistance',   0.17, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'PlotResults',       true, @(x)islogical(x)||ismember(x,[0,1]));
addParameter(p,'PeakOffset',  0, @(x)isnumeric(x)&&isscalar(x)&&mod(x,1)==0);
addParameter(p,'ValleyOffset',0, @(x)isnumeric(x)&&isscalar(x)&&mod(x,1)==0);
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

    % --- 2) Apply offsets safely ---
    locs_peaks   = locs_peaks   + params.PeakOffset;
    locs_valleys = locs_valleys + params.ValleyOffset;

    % Clip to valid range [1, N]
    locs_peaks(locs_peaks < 1 | locs_peaks > N) = [];
    locs_valleys(locs_valleys < 1 | locs_valleys > N) = [];

    % --- 3) Build predicted mask ---
    mask = peakValleyMask(locs_peaks, locs_valleys, N);

    % --- 4) Align from first peak to end ---
    if isempty(locs_peaks)
        warning('No peaks found in segment %d.', i);
        continue;
    end
    startPK = locs_peaks(1);
    mask_pred = mask(startPK:end);
    mask_meas = compressor_chunk(startPK:end);
    mask_meas = mask_meas(1:numel(mask_pred));

    % --- 5) Compare predicted vs measured ---
    stats = compareOnTime(mask_pred, mask_meas);

    % --- 6) Save results ---
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

    % --- 7) Plot results if requested ---
    if params.PlotResults
        figure('Name', sprintf('Peaks/Valleys Segment %d', i))
        plot(t_chunk/3600, x_chunk, 'b-'); hold on
        plot(t_chunk(locs_peaks)/3600, pks(1:numel(locs_peaks)), 'ro', 'MarkerFaceColor', 'r');
        plot(t_chunk(locs_valleys)/3600, valleys(1:numel(locs_valleys)), 'go', 'MarkerFaceColor', 'g');
        xlabel('Time (h)'), ylabel('Amplitude')
        title(sprintf('Peak Detection — Segment %d', i))
        grid on, hold off

        figure('Name', sprintf('Compressor Compare Segment %d', i))
        plot(mask_pred, '-r', 'LineWidth', 2); hold on
        plot(mask_meas, '-b'); hold off
        xlabel('Sample index'), ylabel('Compressor ON/OFF')
        legend('Predicted','Measured')
        title(sprintf('Compressor ON/OFF Comparison — Segment %d', i))
        grid on
    end

    % --- 8) Print stats summary ---
    fprintf('--- Segment %d ---\n', i)
    fprintf('Measured ON-time   = %.1f h\n', stats.timeOnMeas/3600)
    fprintf('Predicted ON-time  = %.1f h\n', stats.timeOnPred/3600)
    fprintf('Difference         = %+0.1f s (%.1f%%)\n\n', ...
            stats.diffTime, stats.pctError)
end
end
