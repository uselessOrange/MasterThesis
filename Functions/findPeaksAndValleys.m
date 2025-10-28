function [pks, locs_peaks, valleys, locs_valleys] = findPeaksAndValleys(data, varargin)
% findPeaksAndValleys   Find peaks and valleys in a 1-D signal.
%
% [pks, locsPeaks, vlys, locsVlys] = findPeaksAndValleys(data)
%   uses the defaults:
%     MinPeakProminence = 0.8
%     MinPeakDistance   = 0.17
%
% [pks, locsPeaks, vlys, locsVlys] = findPeaksAndValleys(data, ...)
%   you may pass ANY of the Name–Value pairs that MATLAB’s findpeaks()
%   supports, e.g. 'MinPeakHeight', H, 'Threshold', T, 'MaxPeakWidth',W, ...
%
% Inputs:
%   data                – 1×N or N×1 vector to search
%   … any Name–Value pairs recognized by findpeaks()
%
% Outputs:
%   pks       – amplitudes of the positive peaks
%   locsPeaks – indices (or X-values) of those peaks in data
%   valleys   – amplitudes of the valleys (in data)
%   locsVlys  – indices (or X-values) of those valleys

  % your defaults
  defaultProm = 0.8;
  defaultDist = 0.17;

  % if user did not supply their own Prominence, inject default
  if ~any(strncmpi('MinPeakProminence', varargin, numel('MinPeakProminence')))
    varargin = [{'MinPeakProminence', defaultProm}, varargin];
  end

  % if user did not supply their own Distance, inject default
  if ~any(strncmpi('MinPeakDistance', varargin, numel('MinPeakDistance')))
    varargin = [{'MinPeakDistance', defaultDist}, varargin];
  end

  % --- now pass ALL name–value pairs on to findpeaks() ---
  [pks,      locs_peaks]   = findpeaks(data,    varargin{:});
  [negVlys, locs_valleys]  = findpeaks(-data,   varargin{:});
  valleys = -negVlys;

end