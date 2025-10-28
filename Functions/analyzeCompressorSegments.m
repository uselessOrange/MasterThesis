function results = analyzeCompressorSegments(t, x, compressor, usableDataIndex, varargin)
% analyzeCompressorSegments
%   Performs peak/valley analysis (via findPeaksAndValleys), compressor ON/OFF
%   mask prediction, and comparison against measured data for multiple segments.
%
%   You may pass ANY of the Name–Value pair arguments that findpeaks()
%   supports (MinPeakHeight, Threshold, MaxPeakWidth, …), *plus* these extras:
%     'PlotResults'  (true/false, default true)
%     'PeakOffset'   (integer,    default 0)
%     'ValleyOffset' (integer,    default 0)
%
% Example:
%   results = analyzeCompressorSegments(t, x, compressor, usableIndex, ...
%                  'PeakOffset',5, 'ValleyOffset',-3, ...
%                  'MinPeakProminence', 1.2, 'MinPeakDistance',0.2);

  %% --- 1) Parse own parameters, keep the rest unmatched ---
  p = inputParser;
  p.KeepUnmatched = true;
  addParameter(p, 'PlotResults',   true,  @(v) islogical(v) || ismember(v,[0,1]));
  addParameter(p, 'PeakOffset',    0,     @(v) isnumeric(v)&&isscalar(v)&&mod(v,1)==0);
  addParameter(p, 'ValleyOffset',  0,     @(v) isnumeric(v)&&isscalar(v)&&mod(v,1)==0);
  parse(p, varargin{:});
  params = p.Results;

  % build a cell array of everything else to forward to findPeaksAndValleys:
  fn = fieldnames(p.Unmatched);
  forwardArgs = cell(1,2*numel(fn));
  for k = 1:numel(fn)
    forwardArgs{2*k-1} = fn{k};
    forwardArgs{2*k  } = p.Unmatched.(fn{k});
  end

  %% --- 2) Initialize output struct ---
  nSets = numel(usableDataIndex);
  results = repmat(struct( ...
     'pks',[], 'locs_peaks',[], 'valleys',[], 'locs_valleys',[], ...
     'mask_pred',[], 'mask_meas',[], 'stats',[], ...
     't_chunk',[], 'x_chunk',[], 'compressor_chunk',[] ), 1, nSets);

  %% --- 3) Loop over segments ---
  for i = 1:nSets
    idx = usableDataIndex{i};
    t_chunk          = t(idx);
    x_chunk          = x(idx);
    compressor_chunk = compressor(idx);
    N = numel(idx);

    % 3.1) call your wrapper, passing along *all* findpeaks name–value pairs
    [pks, locs_peaks, valleys, locs_valleys] = ...
       findPeaksAndValleys(x_chunk, forwardArgs{:});

    % 3.2) apply the integer offsets
    locs_peaks   = locs_peaks   + params.PeakOffset;
    locs_valleys = locs_valleys + params.ValleyOffset;

    % clip out-of-range indices
    validP = locs_peaks   >=1 & locs_peaks   <= N;
    validV = locs_valleys >=1 & locs_valleys <= N;
    locs_peaks   = locs_peaks(validP);
    pks          = pks(validP);
    locs_valleys = locs_valleys(validV);
    valleys      = valleys(validV);

    % 3.3) build predicted ON/OFF mask
    mask = peakValleyMask(locs_peaks, locs_valleys, N);

    % 3.4) align from first peak to end
    if isempty(locs_peaks)
      warning('No peaks found in segment %d.', i);
      continue
    end
    startPK    = locs_peaks(1);
    mask_pred  = mask(startPK:end);
    mask_meas  = compressor_chunk(startPK:end);
    mask_meas  = mask_meas(1:numel(mask_pred));

    % 3.5) compare
    stats = compareOnTime(mask_pred, mask_meas);

    % 3.6) store into results
    results(i).pks              = pks;
    results(i).locs_peaks       = locs_peaks;
    results(i).valleys          = valleys;
    results(i).locs_valleys     = locs_valleys;
    results(i).mask_pred        = mask_pred;
    results(i).mask_meas        = mask_meas;
    results(i).stats            = stats;
    results(i).t_chunk          = t_chunk;
    results(i).x_chunk          = x_chunk;
    results(i).compressor_chunk = compressor_chunk;

    % 3.7) optional plotting
    if params.PlotResults
      figure('Name', sprintf('Peaks/Valleys Segment %d', i))
      plot(t_chunk/3600, x_chunk, 'b-'); hold on
      plot(t_chunk(locs_peaks)/3600, pks,      'ro','MarkerFaceColor','r');
      plot(t_chunk(locs_valleys)/3600, valleys,'go','MarkerFaceColor','g');
      xlabel('Time (h)'), ylabel('Amplitude')
      title(sprintf('Segment %d: Detected Peaks & Valleys', i))
      grid on, hold off

      figure('Name', sprintf('Compressor Compare Segment %d', i))
      plot(mask_pred, '-r','LineWidth',2); hold on
      plot(mask_meas,  '-b'); 
      xlabel('Sample index'), ylabel('ON/OFF')
      legend('Predicted','Measured')
      title(sprintf('Segment %d: ON/OFF Comparison', i))
      grid on
    end

    % 3.8) print stats
    fprintf('--- Segment %d ---\n', i)
    fprintf('Measured ON-time   = %.1f h\n', stats.timeOnMeas/3600)
    fprintf('Predicted ON-time  = %.1f h\n', stats.timeOnPred/3600)
    fprintf('Difference         = %+0.1f s (%.1f%%)\n\n', ...
            stats.diffTime, stats.pctError)
  end
  
end