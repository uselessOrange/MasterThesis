function mask = peakValleyMask(locs_peaks, locs_valleys, N)
% peakValleyMask  returns a length‐N vector that is 1 between each peak
%                and the next valley, and 0 between each valley and the
%                next peak.
%
%  Inputs:
%    locs_peaks   = vector of sample‐indices of your peaks
%    locs_valleys = vector of sample‐indices of your valleys
%    N            = desired length of output mask
%
%  Output:
%    mask = 1×N vector of 0/1

  % Pre-allocate
  mask = zeros(1,N);  

  % Build [index, value] array: value=1 for peaks, 0 for valleys
  P = [locs_peaks(:),   ones(numel(locs_peaks),1)];
  V = [locs_valleys(:), zeros(numel(locs_valleys),1)];
  events = sortrows([P; V],1);

  % Find first peak in the sorted events
  firstPeakIdx = find(events(:,2)==1,1,'first');
  if isempty(firstPeakIdx)
    disp("no peaks at all => all zeros")
    return
  end

  % Now walk through events from that peak to the end
  for k = firstPeakIdx : size(events,1)-1
    thisPos   = events(k,1);
    thisValue = events(k,2);        % 1 if peak, 0 if valley
    nextPos   = events(k+1,1);
    % Fill mask(thisPos : nextPos) with thisValue
    % (clamp nextPos to N so we don’t overrun)
    j1 = max(1,min(N,thisPos));
    j2 = max(1,min(N,nextPos));
    mask(j1:j2) = thisValue;
  end

  % (Optionally) fill beyond the last event with the last event’s value:
  lastPos   = events(end,1);
  lastValue = events(end,2);
  if lastPos < N
    mask(lastPos+1:end) = lastValue;
  end
end