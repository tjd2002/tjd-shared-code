function c = continterp(c,varargin)
% CONTINTERP resample and interpolate cont struct data to match timestamps
%
% param/value pair args:
%
%  'timewin', tstart/tend of resulting cdat (default same as input)
%  'nsamps'/'samplerate': specify interp points
%  'method', interp1 method ({'cubic'}, 'spline', 'linear', 'nearest' etc)
%            (NB: spline will propagate NaNs through all data)
%  'resampbeforeinterp', should we resample near new samplerate before
%  interpolating? (avoids aliasing, but increases edge effects for
%  signals with lots of NaN/Inf)
%
% todo:
%   -performance: bracket timewin before resample&interp1
  
  a = struct(...
      'timewin',[],...
      'nsamps', [],...
      'samplerate',[],...
      'resampbeforeinterp', true,...
      'method', 'cubic');
  
  a = parseArgsLite(varargin,a);
  
  if isempty(a.timewin),
    if size(c.data,1) > 1e6;
      warning('no timewin provided interpolating entire cont struct');
    end

    a.timewin = [c.tstart c.tend];
    
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
    return;
  end
  
  %%% [bracket requested timerange]
  % pad with 2 samples, but don't extend past tstart/tend
  timewin = a.timewin + [-2 2]./(min([c.samplerate a.samplerate]));
  timewin = [max(timewin(1), c.tstart) min(timewin(2), c.tend)];

  c = contwin(c, timewin);
  
  %%% initial downsample of high sampling rate data to avoid aliasing
  if ~strcmp(a.method, 'nearest') 
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
  for k = nchans:-1:1 % reverse order so that first loop preallocates
    % blech: interp1 returns data in a row for vector inputs, 
    % last argument says use NaN for extrapolated values
    newdata(:,k) = interp1(x,c.data(:,k), xi, a.method, NaN);
  end

  c.data = newdata;
  
  c.tstart = xi(1);
  c.tend = xi(end);
  c.samplerate = samplerate_effective;
  c = contdatarange(c);

  % hard to say where the bad samples ended up
  c.nbad_start = NaN;
  c.nbad_end = NaN;