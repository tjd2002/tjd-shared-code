function c_out = continterp(c,varargin)
% CONTINTERP resample and interpolate cont data at specified sample times
%
%  c_out = continterp(c,[name/value pairs])
%
% The sample times at which to interpolate the data are specified by
% providing a start/end time, and either a total number of samples, or a
% new sampling rate.
%
% Inputs: (* = required)
%  *c - cont struct to be interpolated
%  'timewin' - start/end time of resulting cdat (default same as input)
%  *'nsamps'/'samplerate' - number of samples or sampling rate (in Hz) for c_out.
%  'method', interpolation method ({'pchip'}, 'spline', 'linear', 'nearest' etc)
%      (NB: 'spline' not recommended; will propagate NaNs through all 
%      data, result in inaccurate nbad_start/nbad_end estimates)
%
%  Uncommon optional inputs:
%  'resampbeforeinterp', should we resample near new samplerate before
%      interpolating? This is important to avoid aliasing, but when a signal
%      contains short stretches of data surrounded by NaN/Inf values, it
%      can worsen the edge effects caused by interpolation. (Defaults to
%      'true', except that for 'nearest' interpolation, it is always false).
%  'extrapval': argument to interp1. In unusual cases, can end up needing to extrapolate
%      first/last sample. Default is 'Nan', [] means to omit the extrapval 
%      argument, allowing extrapolation when  using 'method' pchip or spline.
%      Including this argument suppresses the warning about extrapolation.
%
%  Outputs:
%  c_out - cont struct with new timebase
  
% Tom Davidson <tjd@stanford.edu> 2003-2010

  % data integrity check
  contcheck(c);

  a = struct(...
      'timewin',[],...
      'nsamps', [],...
      'samplerate',[],...
      'resampbeforeinterp', true,...
      'extrapval', [],...
      'method', 'pchip');
  
  a = parseArgsLite(varargin,a);
  
  if isempty(a.timewin),
% $$$     if size(c.data,1) > 1e6;
% $$$       warning('no timewin provided interpolating entire cont struct');
% $$$     end
    a.timewin = [c.tstart c.tend];
  end
  
  if isempty(a.extrapval),
      extrapwarn = true;
      a.extrapval = NaN;
  else
      extrapwarn = false;
  end
  
  if strcmp(lower(a.method), 'cubic')
    %'cubic' method is deprecated in R2014a, sub in equivalent 'pchip'
    a.method = 'pchip';
  end
    
  if sum([~isempty(a.nsamps) ~isempty(a.samplerate)]) ~= 1,
    error('exactly one of nsamps or samplerate must be provided');
  end

  if ~isempty(a.nsamps)
    samplerate_effective = (a.nsamps-1)./diff(a.timewin);
  else
    samplerate_effective = a.samplerate;
  end
  
  %%% return if no interpolation requested
  if ((~isempty(a.nsamps) && a.nsamps == size(c.data,1)) || ...
      (~isempty(a.samplerate) && a.samplerate == c.samplerate)) ...
        && ...
        a.timewin(1) == c.tstart &&...
        a.timewin(2) == c.tend
    disp('No interpolation needed.');
    c_out = c;
    return;
  end
  
  if a.timewin(1)<c.tstart || a.timewin(2)>c.tend && extrapwarn,
      warning('Requested timewin [%g %g] outside input cdat range [%g %g], data will be padded with ''extrapval''',...
          a.timewin(1), a.timewin(2), c.tstart, c.tend);
  end
  
  %%% [bracket requested timerange for cropping]
  % pad with 2 samples (larger of pre/post samples),
  % but don't extend past tstart/tend
  cropwin = a.timewin + [-2 2]./(min([c.samplerate a.samplerate]));
  cropwin = [max(cropwin(1), c.tstart) min(cropwin(2), c.tend)];

  % crop data early to save filtering time
  c = contwin(c, cropwin);
  
  %%% initial downsample of high sampling rate data to avoid aliasing
  if strcmp(a.method, 'nearest') 
    warning([mfilename ':NearestMethod'], ...
            ['''nearest'' method implies no resampling before interpolation. ' ...
             'Beware aliasing']);
    a.resampbeforeinterp = false;
  end
  if samplerate_effective < c.samplerate  && a.resampbeforeinterp
    disp('Resampling before interpolation...');
    c = contresamp(c, ...
                   'resample', samplerate_effective./c.samplerate,...
                   'tol', 0.05); % use wide tolerance since we're going to
                                 % interp below, anyway
  
  else
    disp('no resamp before interp');
  end
  
  %%% interpolate!

  % generate timestamps of existing data
  x = linspace(c.tstart,c.tend, size(c.data,1));
  
  % generate timestamps of requested samples
  if ~isempty(a.nsamps),
    xi = linspace(a.timewin(1), a.timewin(2), a.nsamps);
  else
    xi = a.timewin(1):1/a.samplerate:a.timewin(2);
  end
  
  disp('interpolating...');

  nchans = size(c.data,2);

  % This causes these values to get interpolated to Nan (otherwise interp1
  % *will* interpolate over internal NaNs; the 'Nan' extrapval argument
  % to interp1 only applies to points outside x)
  c.data(isnan(c.data)) = -Inf;
  
  % do it channel-by-channel to save memory
  for k = nchans:-1:1 % reverse loop so that first iteration preallocates
    % blech: interp1 returns data in a row for vector inputs, 
    % last argument says use NaN for extrapolated values
    if ~isempty(a.extrapval)
      newdata(:,k) = interp1(x,c.data(:,k), xi, a.method, a.extrapval);
    else
      newdata(:,k) = interp1(x,c.data(:,k), xi, a.method);
    end
  end
   
  c_out = c; % initialize with old cont struct
  c_out.data = newdata;
  
  c_out.tstart = xi(1);
  c_out.tend = xi(end);
  c_out.samplerate = samplerate_effective;

  % Avoid attempting to crop past end of data in case where there was no
  % padding data (e.g. when interping to end of data)
  if a.timewin(2) > c_out.tend,
      a.timewin(2) = c_out.tend;
  end
  
  % recrop data more tightly using a.timewin (remove any timewin padding)
  c_out = contwin(c_out, a.timewin, 'samps_nearest');

  % calculate nbad_start/_end
  
  % find *time* of good/bad samples (after any possible filtering in this 
  % function).
  timewin_good(1) = c.tstart + (c.nbad_start-1)./c.samplerate;
  timewin_good(2) = c.tend - (c.nbad_end-1)./c.samplerate;
  
  % Excldue extra one sample (in underlying interp'd signal) to deal with 
  % effects of interpolation at edge of bad data. NB: Not totally 
  % conservative; e.g. pchip effects could last >1 sample but seems OK.
  timewin_good = timewin_good + ([1 -1] ./ c.samplerate);
  
  % Use these *times* to find nbad_start/end in output (if any)
  % (Bracket timewin with floor/ceil because we are reporting the number
  % of *bad* samples)
  c_out.nbad_start = max(0, floor( c_out.samplerate .* (timewin_good(1) - c_out.tstart)));
  c_out.nbad_end   = max(0, ceil( c_out.samplerate .* (c_out.tend - timewin_good(2))));

  
  % may not be necessary?
  c_out = contdatarange(c_out);
  
  % data integrity check
  contcheck(c_out);