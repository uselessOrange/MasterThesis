function [stats] = compareOnTime(mask_pred, mask_meas, dt)
% compareOnTime   compare total ON‐time of prediction vs. measurement
%
%   Inputs:
%     mask_pred = 1×N logical or numeric 0/1 vector (predicted compressor on/off)
%     mask_meas = 1×N logical or numeric 0/1 vector (measured compressor on/off)
%     dt        = (optional) sampling interval in seconds (default = 1)
%
%   Output:
%     stats is a struct with fields
%       .nOnPred    = total ON samples predicted
%       .nOnMeas    = total ON samples measured
%       .timeOnPred = total ON time predicted (seconds)
%       .timeOnMeas = total ON time measured (seconds)
%       .diffSamp   = nOnPred − nOnMeas (samples)
%       .diffTime   = timeOnPred − timeOnMeas (seconds)
%       .absDiff    = |.diffTime| (seconds)
%       .pctError   = 100*(timeOnPred − timeOnMeas)/timeOnMeas (%)

  if nargin<3
    dt = 1;
  end

  if numel(mask_pred)~=numel(mask_meas)
    error('Prediction and measurement must be the same length.');
  end

  % force logical
  mask_pred = mask_pred(:)~=0;
  mask_meas = mask_meas(:)~=0;

  % Counts
  stats.nOnPred  = sum(mask_pred);
  stats.nOnMeas  = sum(mask_meas);

  % Convert to time
  stats.timeOnPred = stats.nOnPred * dt;
  stats.timeOnMeas = stats.nOnMeas * dt;

  % Differences
  stats.diffSamp = stats.nOnPred - stats.nOnMeas;
  stats.diffTime = stats.timeOnPred - stats.timeOnMeas;
  stats.absDiff  = abs(stats.diffTime);

  % Percent error relative to measured ON‐time
  if stats.timeOnMeas>0
    stats.pctError = 100 * stats.diffTime / stats.timeOnMeas;
  else
    stats.pctError = NaN;
  end
end