function results = analyzeCompressorMaskFromGrayBox(t, x_mask, compressor, usableDataIndex, varargin)
% analyzeCompressorMaskSegments
%   Performs ON/OFF mask comparison between a predicted compressor mask (x_mask)
%   and the measured compressor signal for multiple segments.
%
%   This function is similar to analyzeCompressorSegments, except that the
%   compressor ON/OFF mask is provided directly (x_mask), so no peak/valley
%   detection is performed.
%
% Example:
%   results = analyzeCompressorMaskSegments(t, predMask, compressor, usableIndex, ...
%                  'PlotResults', true);

  %% --- 1) Parse own parameters ---
  p = inputParser;
  addParameter(p, 'PlotResults', true, @(v) islogical(v) || ismember(v,[0,1]));
  parse(p, varargin{:});
  params = p.Results;

  %% --- 2) Initialize output struct ---
  nSets = numel(usableDataIndex);
  results = repmat(struct( ...
     'mask_pred',[], 'mask_meas',[], 'stats',[], ...
     't_chunk',[], 'mask_chunk',[], 'compressor_chunk',[] ), 1, nSets);

  %% --- 3) Loop over segments ---
  for i = 1:nSets
    idx = usableDataIndex{i};
    t_chunk          = t(idx);
    mask_chunk       = x_mask(idx);
    compressor_chunk = compressor(idx);

    % Ensure same length
    N = min(numel(mask_chunk), numel(compressor_chunk));
    mask_pred = mask_chunk(1:N);
    mask_meas = compressor_chunk(1:N);

    % Compare ON-time and other metrics
    stats = compareOnTime(mask_pred, mask_meas);

    % Store results
    results(i).mask_pred        = mask_pred;
    results(i).mask_meas        = mask_meas;
    results(i).stats            = stats;
    results(i).t_chunk          = t_chunk;
    results(i).mask_chunk       = mask_chunk;
    results(i).compressor_chunk = compressor_chunk;

    % Optional plotting
    if params.PlotResults
      figure('Name', sprintf('Compressor Mask Comparison Segment %d', i))
      plot(mask_pred, '-r','LineWidth',2); hold on
      plot(mask_meas,  '-b');
      xlabel('Sample index'), ylabel('ON/OFF')
      legend('Predicted','Measured')
      title(sprintf('Segment %d: Compressor Mask Comparison', i))
      grid on
    end

    % Print stats
    fprintf('--- Segment %d ---\n', i)
    fprintf('Measured ON-time   = %.1f h\n', stats.timeOnMeas/3600)
    fprintf('Predicted ON-time  = %.1f h\n', stats.timeOnPred/3600)
    fprintf('Difference         = %+0.1f s (%.1f%%)\n\n', ...
            stats.diffTime, stats.pctError)
  end
end
