function segments = computeCompressorStats( ...
    t_day, x_day_interp, compressor_day, usableDataIndex, pkOpts, pvOpts )
% computeCompressorStats
%   For each block in usableDataIndex:
%     – extract time, signal, compressor
%     – detect peaks & valleys
%     – build on/off mask
%     – align to first peak, compare to measured
%   No plotting. Returns a struct array 'segments' with all data.
%
% USAGE:
%   segments = computeCompressorStats( ...
%       t_day, x_day_interp, compressor_day, usableDataIndex );
%   or with custom findpeaks options
%   segments = computeCompressorStats( ...
%       t_day, x_day_interp, compressor_day, usableDataIndex, ...
%       struct('MinPeakProminence',0.1), ...
%       struct('MinPeakProminence',0.1) );
%
% INPUTS:
%   t_day            – time vector (seconds)
%   x_day_interp     – interpolated signal
%   compressor_day   – measured on/off (0/1)
%   usableDataIndex  – 1×M cell array of index vectors
%   pkOpts (opt)     – struct of options passed to findpeaks for peaks
%   pvOpts (opt)     – struct of options passed to findpeaks for valleys
%
% OUTPUT:
%   segments         – 1×M struct array with fields:
%       .idx          – original indices
%       .t_hr         – t_day(idx)/3600
%       .x            – x_day_interp(idx)
%       .comp         – compressor_day(idx)
%       .pks, .locsP  – peak values & locations
%       .valleys, .locsV – valley values & locations
%       .mask         – on/off mask full length N
%       .startPt      – index of first peak
%       .mask_pred    – mask(startPt:end)
%       .mask_meas    – comp(startPt:end) (trimmed to same length)
%       .stats        – compareOnTime result
%
if nargin<5, pkOpts = struct(); end
if nargin<6, pvOpts = struct(); end

nSeg = numel(usableDataIndex);
segments = repmat(struct( ...
  'idx',[], 't_hr',[], 'x',[], 'comp',[], ...
  'pks',[], 'locsP',[], 'valleys',[], 'locsV',[], ...
  'mask',[], 'startPt',[], 'mask_pred',[], 'mask_meas',[], ...
  'stats',[] ), 1, nSeg);

for k = 1:nSeg
  idx = usableDataIndex{k};
  t_k = t_day(idx)/3600;
  x_k = x_day_interp(idx);
  c_k = compressor_day(idx);
  N   = numel(idx);

  % find peaks
  [pks, locsP] = findpeaks(x_k, pkOpts);
  % find valleys
  [nv, locsV] = findpeaks(-x_k, pvOpts);
  valleys = -nv;

  % build on/off mask
  mask = peakValleyMask(locsP, locsV, N);

  % align from first peak
  startPt    = locsP(1);
  mp         = mask(startPt:end);
  mm         = c_k(startPt:end);
  L          = min(numel(mp), numel(mm));
  mp = mp(1:L);
  mm = mm(1:L);

  % compare
  stats = compareOnTime(mp, mm);

  % store
  segments(k).idx       = idx;
  segments(k).t_hr      = t_k;
  segments(k).x         = x_k;
  segments(k).comp      = c_k;
  segments(k).pks       = pks;
  segments(k).locsP     = locsP;
  segments(k).valleys   = valleys;
  segments(k).locsV     = locsV;
  segments(k).mask      = mask;
  segments(k).startPt   = startPt;
  segments(k).mask_pred = mp;
  segments(k).mask_meas = mm;
  segments(k).stats     = stats;
end
end