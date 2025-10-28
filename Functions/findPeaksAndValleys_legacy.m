function [pks, locs_peaks, valleys, locs_valleys] = findPeaksAndValleys(data, varargin)
% findPeaksAndValleys   Find peaks and valleys in a 1-D signal.
%
%   [pks, locsPeaks, vlys, locsVlys] = findPeaksAndValleys(data)
%   uses the defaults:
%       MinPeakProminence = 0.8
%       MinPeakDistance   = 0.17
%
%   [pks, locsPeaks, vlys, locsVlys] = findPeaksAndValleys(data,...
%                        'MinPeakProminence',P,'MinPeakDistance',D)
%   allows you to override those defaults.
%
% Inputs:
%   data                – 1×N or N×1 vector to search
%   Name-Value Pairs:
%     'MinPeakProminence' (scalar, default 0.8)
%     'MinPeakDistance'   (scalar, default 0.17)
%
% Outputs:
%   pks       – amplitudes of the positive peaks
%   locsPeaks – indices of those peaks in data
%   vlys      – amplitudes of the valleys (negative peaks)
%   locsVlys  – indices of those valleys in data

  % parse inputs
  p = inputParser;
  addRequired (p, 'data',            @(x) isnumeric(x) && isvector(x));
  addParameter(p, 'MinPeakProminence', 0.8, @(x) isnumeric(x) && isscalar(x));
  addParameter(p, 'MinPeakDistance',   0.17, @(x) isnumeric(x) && isscalar(x));
  parse(p, data, varargin{:});

  mp = p.Results.MinPeakProminence;
  md = p.Results.MinPeakDistance;

  % find positive peaks
  [pks, locs_peaks] = findpeaks(data, ...
                               'MinPeakProminence', mp, ...
                               'MinPeakDistance',   md);

  % find valleys by inverting the signal
  [vlysNeg, locs_valleys] = findpeaks(-data, ...
                               'MinPeakProminence', mp, ...
                               'MinPeakDistance',   md);
  valleys = -vlysNeg;

end