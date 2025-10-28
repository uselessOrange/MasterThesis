function [N_badResults, tot_timeUnderConsideration] = ...
    evaluateDataLengthEffects(indexContainer, t, x, compressor, fs, ...
                              runs, TimeScale, varargin)
% evaluateDataLengthEffects  
%   Evaluate # bad predictions vs minimal data‐length.
%
% [N_badResults, tot_timeUnderConsideration] = ...
%     evaluateDataLengthEffects(indexContainer, t, x, compressor, fs, ...
%                               runs, TimeScale, ...)
%
% Optional Name–Value pairs passed *through* to analyzeCompressorSegments
% (and ultimately to findPeaksAndValleys/findpeaks):
%   'PeakOffset'    integer offset (default: 0)
%   'ValleyOffset'  integer offset (default: 0)
% plus *any* of the findpeaks name–value pairs:
%   'MinPeakHeight', 'Threshold', 'MinPeakProminence', …
%
% The call to analyzeCompressorSegments is always made with PlotResults=false,
% but if you supply 'PlotResults',true in varargin it will override that.

  %% 1) Parse only PeakOffset/ValleyOffset, keep all others unmatched
  p = inputParser;
  p.KeepUnmatched = true;
  addParameter(p, 'PeakOffset',   0, @(v) isnumeric(v)&&isscalar(v)&&mod(v,1)==0);
  addParameter(p, 'ValleyOffset', 0, @(v) isnumeric(v)&&isscalar(v)&&mod(v,1)==0);
  parse(p, varargin{:});
  PeakOffset   = p.Results.PeakOffset;
  ValleyOffset = p.Results.ValleyOffset;

  % build a cell array of the unmatched Name–Value pairs
  fn = fieldnames(p.Unmatched);
  forwardNV = cell(1,2*numel(fn));
  for k = 1:numel(fn)
    forwardNV{2*k-1} = fn{k};
    forwardNV{2*k  } = p.Unmatched.(fn{k});
  end

  %% 2) Pre‐allocate outputs
  N_badResults               = nan(1,runs);
  tot_timeUnderConsideration = nan(1,runs);

  %% 3) Loop over minimal‐length steps
  for Min_dataLength = 1:runs
    % build list of segments longer than threshold
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

    % if none left, record zeros
    if isempty(usableDataIndex)
      N_badResults(Min_dataLength) = 0;
      tot_timeUnderConsideration(Min_dataLength) = 0;
      continue
    end

    %% 4) Call analyzer, forcing PlotResults=false
    results = analyzeCompressorSegments( ...
       t, x, compressor, usableDataIndex, ...
       'PeakOffset',   PeakOffset, ...
       'ValleyOffset', ValleyOffset, ...
       'PlotResults',  false, ...
       forwardNV{:} );

    % collect percent‐error from each segment
    pctError = [results.stats];
    pctError = [pctError.pctError];

    % count “bad” (>10%)
    idxBad = abs(pctError) > 10;
    N_badResults(Min_dataLength) = sum(idxBad);
    tot_timeUnderConsideration(Min_dataLength) = totalLen;
  end

  %% 5) Finally, plot summary
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