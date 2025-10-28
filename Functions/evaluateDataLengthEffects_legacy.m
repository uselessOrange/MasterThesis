function [N_badResults, tot_timeUnderConsideration] = ...
    evaluateDataLengthEffects(indexContainer, t, x, compressor, fs, ...
                              runs, TimeScale, varargin)
%EVALUATEDATALENGTHEFFECTS  Evaluate # bad predictions vs minimal data length.
%
%   [N_badResults, tot_timeUnderConsideration] =
%     evaluateDataLengthEffects(indexContainer, t, x, compressor, fs, ...
%                               runs, TimeScale)
%   uses default PeakOffset = 0, ValleyOffset = 0.
%
%   [N_badResults, tot_timeUnderConsideration] =
%     evaluateDataLengthEffects(..., 'PeakOffset', PO, 'ValleyOffset', VO)
%   overrides those defaults.
%
%   Inputs:
%     indexContainer   – cell array of index vectors
%     t, x, compressor – data passed to analyzeCompressorSegments
%     fs               – sampling frequency (Hz)
%     runs             – # of different minimal‐length steps to try
%     TimeScale        – scale from step index to hours
%
%   Optional name-value pairs:
%     'PeakOffset'     – offset passed to analyzeCompressorSegments
%                        (default: 0)
%     'ValleyOffset'   – offset passed to analyzeCompressorSegments
%                        (default: 0)
%
%   Outputs:
%     N_badResults             – 1×runs vector of counts of bad preds (>10%)
%     tot_timeUnderConsideration – 1×runs vector of total samples used

  %— parse optional inputs
  p = inputParser;
  addParameter(p,'PeakOffset',0,@(x) isnumeric(x)&&isscalar(x));
  addParameter(p,'ValleyOffset',0,@(x) isnumeric(x)&&isscalar(x));
  parse(p,varargin{:});
  PeakOffset   = p.Results.PeakOffset;
  ValleyOffset = p.Results.ValleyOffset;

  %— preallocate outputs
  N_badResults               = nan(1,runs);
  tot_timeUnderConsideration = nan(1,runs);

  %— loop over minimal‐data‐length steps
  for Min_dataLength = 1:runs
    % build the list of usable indices
    usableDataIndex = {};
    totalLen = 0;
    set = 1;
    threshold = 60*60*Min_dataLength*TimeScale*fs;
    for i = 1:numel(indexContainer)
      L = numel(indexContainer{i});
      if L > threshold
        usableDataIndex{set} = indexContainer{i};
        totalLen = totalLen + L;
        set = set + 1;
      end
    end

    % if none survive, record zeros and continue
    if isempty(usableDataIndex)
      N_badResults(Min_dataLength) = 0;
      tot_timeUnderConsideration(Min_dataLength) = 0;
      continue
    end

    % call your analyzer
    results = analyzeCompressorSegments( ...
       t, x, compressor, usableDataIndex, ...
       'PeakOffset',  PeakOffset, ...
       'ValleyOffset',ValleyOffset, ...
       'MinPeakProminence',0.8, ...
       'MinPeakDistance', 0.17, ...
       'PlotResults', false);

    % collect % errors
    pctError = nan(1,numel(results));
    for k = 1:numel(results)
      pctError(k) = results(k).stats.pctError;
    end

    % count “bad” (>10%)
    idxBad = abs(pctError) > 10;
    N_badResults(Min_dataLength) = sum(idxBad);
    tot_timeUnderConsideration(Min_dataLength) = totalLen;
  end

  %— optionally plot results
  figure;
  stem((1:runs)*TimeScale, N_badResults, 'filled');
  xlabel('Minimal data length [h]');
  ylabel('Number of bad predictions');
  title('Error >10% vs minimal data sampling length');

  figure;
  stem((1:runs)*TimeScale, tot_timeUnderConsideration, 'filled');
  xlabel('Minimal data length [h]');
  ylabel('Total number of samples used');
  title('Samples used vs minimal data sampling length');
end