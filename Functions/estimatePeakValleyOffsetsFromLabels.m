function calib = estimatePeakValleyOffsetsFromLabels(t, x, compressor, usableDataIndex, varargin)
% estimatePeakValleyOffsetsFromLabels
%   Estimate fixed PeakOffset / ValleyOffset from labeled data by comparing
%   detected peaks/valleys against actual compressor state transitions.
%
%   Intended workflow:
%     1) Run this on labeled data (where compressor ON/OFF is known)
%     2) Get optimal fixed offsets
%     3) Reuse those offsets later on unlabeled data in analyzeCompressorSegments
%
%   METHOD
%   ------
%   For each segment:
%     - detect peaks and valleys in x
%     - detect measured ON/OFF transitions in compressor
%     - match each predicted event to the nearest measured transition
%     - compute time/index difference in samples
%   Then aggregate all differences across all segments:
%     PeakOffset   = robust center of peak differences
%     ValleyOffset = robust center of valley differences
%
%   IMPORTANT
%   ---------
%   You must specify how peaks/valleys correspond to compressor transitions:
%
%     Option 1:
%       peaks   -> ON transitions   (0 -> 1)
%       valleys -> OFF transitions  (1 -> 0)
%
%     Option 2:
%       peaks   -> OFF transitions
%       valleys -> ON transitions
%
%   Use the 'PeakMatches' and 'ValleyMatches' parameters.
%
%   INPUTS
%   ------
%   t               : full time vector
%   x               : full signal vector used for peak/valley detection
%   compressor      : full labeled compressor vector (0/1 or logical)
%   usableDataIndex : cell array of segment indices
%
%   NAME-VALUE PARAMETERS
%   ---------------------
%   'PeakMatches'       : 'on' or 'off'   (default 'on')
%   'ValleyMatches'     : 'on' or 'off'   (default 'off')
%   'MaxMatchDistance'  : max allowed event matching distance in samples
%                         (default Inf)
%   'CenterMethod'      : 'median', 'mean', or 'trimmedmean'
%                         (default 'median')
%   'TrimPercent'       : used if CenterMethod='trimmedmean' (default 20)
%   'PlotResults'       : true/false (default false)
%
%   Any additional unmatched name-value pairs are forwarded to
%   findPeaksAndValleys().
%
%   OUTPUT
%   ------
%   calib is a struct with fields:
%     .PeakOffset
%     .ValleyOffset
%     .PeakDiffs
%     .ValleyDiffs
%     .nPeakMatches
%     .nValleyMatches
%     .perSegment
%     .settings
%
%   EXAMPLE
%   -------
%   calib = estimatePeakValleyOffsetsFromLabels( ...
%       t, x, compressor, usableDataIndex, ...
%       'PeakMatches', 'on', ...
%       'ValleyMatches', 'off', ...
%       'MaxMatchDistance', 50, ...
%       'CenterMethod', 'median', ...
%       'MinPeakProminence', 1.2);
%
%   % Then later:
%   % results = analyzeCompressorSegments(...,
%   %     'PeakOffset', calib.PeakOffset, ...
%   %     'ValleyOffset', calib.ValleyOffset);

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;

addParameter(p, 'PeakMatches', 'on', @(s) ischar(s) || isstring(s));
addParameter(p, 'ValleyMatches', 'off', @(s) ischar(s) || isstring(s));
addParameter(p, 'MaxMatchDistance', Inf, @(v) isnumeric(v) && isscalar(v) && v >= 0);
addParameter(p, 'CenterMethod', 'median', @(s) ischar(s) || isstring(s));
addParameter(p, 'TrimPercent', 20, @(v) isnumeric(v) && isscalar(v) && v >= 0 && v < 100);
addParameter(p, 'PlotResults', false, @(v) islogical(v) || ismember(v,[0 1]));

parse(p, varargin{:});
params = p.Results;

peakMatches   = lower(string(params.PeakMatches));
valleyMatches = lower(string(params.ValleyMatches));
centerMethod  = lower(string(params.CenterMethod));

if ~ismember(peakMatches, ["on","off"])
    error('PeakMatches must be ''on'' or ''off''.');
end
if ~ismember(valleyMatches, ["on","off"])
    error('ValleyMatches must be ''on'' or ''off''.');
end
if ~ismember(centerMethod, ["median","mean","trimmedmean"])
    error('CenterMethod must be ''median'', ''mean'', or ''trimmedmean''.');
end

% Forward all unmatched args to findPeaksAndValleys
fn = fieldnames(p.Unmatched);
forwardArgs = cell(1, 2*numel(fn));
for k = 1:numel(fn)
    forwardArgs{2*k-1} = fn{k};
    forwardArgs{2*k}   = p.Unmatched.(fn{k});
end

%% Initialize
nSeg = numel(usableDataIndex);

perSegment = repmat(struct( ...
    'segmentIndex', [], ...
    'locs_peaks', [], ...
    'locs_valleys', [], ...
    'onIdx', [], ...
    'offIdx', [], ...
    'peakDiffs', [], ...
    'valleyDiffs', [], ...
    'peakOffsetSeg', NaN, ...
    'valleyOffsetSeg', NaN), 1, nSeg);

allPeakDiffs = [];
allValleyDiffs = [];

%% Process each segment
for i = 1:nSeg
    idx = usableDataIndex{i};
    x_chunk = x(idx);
    c_chunk = compressor(idx);

    % Detect peaks and valleys
    [~, locs_peaks, ~, locs_valleys] = findPeaksAndValleys(x_chunk, forwardArgs{:});

    % Measured transitions from compressor label
    cbin = c_chunk(:) ~= 0;
    dc = diff(cbin);
    onIdx  = find(dc == 1) + 1;   % 0 -> 1
    offIdx = find(dc == -1) + 1;  % 1 -> 0

    % Match peaks
    if peakMatches == "on"
        peakTarget = onIdx;
    else
        peakTarget = offIdx;
    end

    % Match valleys
    if valleyMatches == "on"
        valleyTarget = onIdx;
    else
        valleyTarget = offIdx;
    end

    peakDiffs = nearestMatchedDiffs(locs_peaks, peakTarget, params.MaxMatchDistance);
    valleyDiffs = nearestMatchedDiffs(locs_valleys, valleyTarget, params.MaxMatchDistance);

    allPeakDiffs = [allPeakDiffs, peakDiffs];
    allValleyDiffs = [allValleyDiffs, valleyDiffs];

    perSegment(i).segmentIndex    = i;
    perSegment(i).locs_peaks      = locs_peaks;
    perSegment(i).locs_valleys    = locs_valleys;
    perSegment(i).onIdx           = onIdx;
    perSegment(i).offIdx          = offIdx;
    perSegment(i).peakDiffs       = peakDiffs;
    perSegment(i).valleyDiffs     = valleyDiffs;
    perSegment(i).peakOffsetSeg   = localCenter(peakDiffs, centerMethod, params.TrimPercent);
    perSegment(i).valleyOffsetSeg = localCenter(valleyDiffs, centerMethod, params.TrimPercent);

    if params.PlotResults
        figure('Name', sprintf('Offset calibration segment %d', i));

        subplot(2,1,1)
        stem(locs_peaks, ones(size(locs_peaks)), 'r', 'filled'); hold on
        stem(locs_valleys, -ones(size(locs_valleys)), 'g', 'filled');
        stem(onIdx, 0.6*ones(size(onIdx)), 'b');
        stem(offIdx, -0.6*ones(size(offIdx)), 'k');
        grid on
        xlabel('Sample index')
        ylabel('Event type')
        title(sprintf('Segment %d events', i))
        legend('Peaks','Valleys','Measured ON','Measured OFF')

        subplot(2,1,2)
        hold on
        if ~isempty(peakDiffs)
            plot(peakDiffs, 'ro-');
        end
        if ~isempty(valleyDiffs)
            plot(valleyDiffs, 'go-');
        end
        yline(0,'k--');
        grid on
        xlabel('Matched event #')
        ylabel('Measured idx - predicted idx (samples)')
        title('Matched event offsets')
        legend('Peak diffs','Valley diffs','zero')
    end
end

%% Global offsets
PeakOffset   = localCenter(allPeakDiffs, centerMethod, params.TrimPercent);
ValleyOffset = localCenter(allValleyDiffs, centerMethod, params.TrimPercent);

%% Output
calib = struct();
calib.PeakOffset     = round(PeakOffset);
calib.ValleyOffset   = round(ValleyOffset);
calib.PeakDiffs      = allPeakDiffs(:)';
calib.ValleyDiffs    = allValleyDiffs(:)';
calib.nPeakMatches   = numel(allPeakDiffs);
calib.nValleyMatches = numel(allValleyDiffs);
calib.perSegment     = perSegment;
calib.settings       = params;

fprintf('Estimated offsets from labeled data:\n');
fprintf('  PeakOffset   = %d samples\n', calib.PeakOffset);
fprintf('  ValleyOffset = %d samples\n', calib.ValleyOffset);
fprintf('  Peak matches   : %d\n', calib.nPeakMatches);
fprintf('  Valley matches : %d\n', calib.nValleyMatches);

end


function diffs = nearestMatchedDiffs(predIdx, measIdx, maxMatchDistance)
% nearestMatchedDiffs
%   For each predicted event index, find nearest measured event index and
%   return difference:
%       diff = measured_index - predicted_index
%
%   Only keep matches within maxMatchDistance.
%
%   Note:
%   This allows repeated use of the same measured event by multiple
%   predicted events. That is simple and often acceptable for calibration.
%   If needed later, this can be changed to one-to-one matching.

predIdx = predIdx(:);
measIdx = measIdx(:);

diffs = [];

if isempty(predIdx) || isempty(measIdx)
    return;
end

for i = 1:numel(predIdx)
    [dmin, j] = min(abs(measIdx - predIdx(i)));
    if ~isempty(j) && dmin <= maxMatchDistance
        diffs(end+1) = measIdx(j) - predIdx(i); %#ok<AGROW>
    end
end
end


function c = localCenter(x, method, trimPercent)
% localCenter
%   Robust center estimate for offset aggregation.

if isempty(x)
    c = NaN;
    return;
end

switch lower(string(method))
    case "median"
        c = median(x, 'omitnan');
    case "mean"
        c = mean(x, 'omitnan');
    case "trimmedmean"
        c = trimmean(x, trimPercent);
    otherwise
        error('Unknown center method.');
end
end