function results = analyzeCompressorMaskFromGrayBox(compare_container, varargin)
% analyzeCompressorMaskFromGrayBox
%   Compares a predicted compressor mask (x_mask) against the measured
%   compressor signal for multiple segments without performing peak/valley
%   detection. This function directly uses the provided ON/OFF mask for
%   analysis and outputs relevant statistics for each segment.
%
%   Inputs:
%       compare_container - Struct containing predicted and measured masks.
%       varargin - Optional parameters (e.g., 'PlotResults', true/false).
%
%   Outputs:
%       results - Struct array containing the predicted mask, measured mask,
%                 statistics.
%
% Example:
%   results = analyzeCompressorMaskFromGrayBox(t, predMask, compressor, usableIndex, ...
%                  'PlotResults', true);

  %% --- 1) Parse own parameters ---
  p = inputParser;
  addParameter(p, 'PlotResults', true, @(v) islogical(v) || ismember(v,[0,1]));
  parse(p, varargin{:});
  params = p.Results;

  %% --- 2) Initialize output struct ---
  nSets = numel(compare_container.mask_meas);
  results = repmat(struct( ...
     'mask_pred',[], 'mask_meas',[], 'stats',[], ...
     'mask_chunk',[], 'compressor_chunk',[] ), 1, nSets);

  %% --- 3) Loop over segments ---
  for i = 1:nSets
    mask_chunk       = compare_container.mask_pred{i};
    compressor_chunk = compare_container.mask_meas{i};

    % Ensure same length
    N = min(numel(mask_chunk), numel(compressor_chunk));
    if numel(mask_chunk) ~= numel(compressor_chunk)
        error('Error: The predicted mask and measured mask must be of the same length.');
    end
    mask_pred = mask_chunk(1:N);
    mask_meas = compressor_chunk(1:N);

    % Compare ON-time and other metrics
    stats = compareOnTime(mask_pred, mask_meas);

    % Store results
    results(i).mask_pred        = mask_pred;
    results(i).mask_meas        = mask_meas;
    results(i).stats            = stats;
    results(i).mask_chunk       = mask_chunk;
    results(i).compressor_chunk = compressor_chunk;

    % Optional plotting
    if params.PlotResults
      figure('Name', sprintf('Compressor Mask Comparison Segment %d', i))
      plot(mask_pred, '-r','LineWidth',2); hold on
      plot(mask_meas,  '-b');hold off
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
