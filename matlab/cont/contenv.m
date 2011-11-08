function c = contenv (c,varargin)
% CONTENV - Compute instantaneous amplitude for data in a cont struct
%
%    cout = contenv(c, [name/value pair args]);
%
% Calculates the instantaneous magnitude (or ENVelope) of a signal, by one of
% several methods (see code for details of algorithms):
%
%  -'peaks' rectifies the signal, and performs linear interpolation across
%   the peaks of the resulting--closest to an envelope.
%  -'hilbert' calculates the magnitude of the analytic signal, which is
%   equal to the amplitude envelope for narrow bandwidth signals
%   (e.g. those which have been bandpass filtered in a narrow
%   range). ('hilbert_complex' returns the complex, analytic signal instead
%   of the envelope)
%  -'rms' calculates the root-mean-square amplitude, smoothed with a
%   sliding boxcar/rectangular filter. 
%
% Inputs: (* means required, -> indicates default value)
%  * c - input cont struct 
%   'method' - method to use to compute envelope.
%       (->'peaks', 'hilbert', 'hilbert_complex', 'rms')
%   'rms_window_t' - time to average if using the 'rms' method (no default)
%
%  Infrequently-used options
%   'envopt' - envelope options structure as created by mkenvopt
%
% Outputs:
%  cout - cont struct with the envelope
%
% Example:
%  cout = contenv(cdat, 'method', 
%
% Tom Davidson <tjd@alum.mit.edu> 2003-2010 

  % data integrity check
  contcheck(c);
  
  a = struct(...
      'envopt',[],...
      'method',[],...
      'rms_window_t',[],...
      'nosuffix', false);
  
  a = parseArgsLite(varargin,a);

  if ~isempty(a.envopt) && ~isempty(a.method),
    error('Can''t provide both ''method'' and ''envopt'' arguments')
  end
  
  if isempty(a.envopt) && isempty(a.method),
    disp('no envelope method requested, using mkenvopt defaults')
    a.envopt = mkenvopt;
  end
  
  if ~isempty(a.method)
    a.envopt = mkenvopt('method', a.method, ...
                        'rms_window_t', a.rms_window_t);
  end
  
  [nsamps nchans] = size(c.data); %#ok
  
  disp('calculating...');

  switch(a.envopt.method)
   case 'hilbert'
    suffix = '_env_hilb';

    if any(isnan(c.data))
      error(['Cannot calculate Hilbert transform on data with NaNs; fix ' ...
             'data or use another method']); %#ok
    end

    for k = 1:nchans,

      % hilbert seems to be buggy with 'single' data; cast to double and back
      dtype = class(c.data);
      
      % get amplitude envelope (magnitude of the analytic signal)
      tmp = abs(hilbert(double(c.data(:,k))));

      % store result in the original data type
      c.data(:,k) = cast(tmp, dtype);

    end
    
   case 'peaks',
    suffix = '_env_pks';

    for k = 1:nchans,
      % rectify the channel
      abschan = abs(c.data(:,k));
      % find local maxima (peaks)
      pks_idx = localmax(abschan);

      if sum(pks_idx>=2)
        % linear interpolation between peaks
        c.data(:,k) = interp1q(find(pks_idx(:,k)), ...
                               abschan(pks_idx), ...
                               (1:nsamps)');
        
      else
        % signal is constant or only has 1 peak--no valid interpolation
        % between peaks
        c.data(:,k) = NaN(size(c.data(:,k)));

      end
    end
    
   case 'rms', 
    if isempty(a.envopt.rms_window_t),
      error('''rms_window_t'' must be provided for ''rms'' method');
    end
    suffix = ['_rms_' num2str(a.envopt.rms_window_t*1000) 'ms'];
    
    % save current chanlabels;
    oldcl = c.chanlabels;
    
    if a.envopt.rms_window_t < (1.5/c.samplerate)
      warning(['rms window will be 1 sample or less, returning original ' ...
               'signal']);
    else
      % do RMS:
      % 1) Square signal
      c.data = c.data.^2;
      % 2) calculate moving mean of squared signal
      c = contfilt(c, 'filtopt', mkfiltopt('name', 'rms_averaging',...
                                           'filttype', 'rectwin',...
                                           'length_t', a.envopt.rms_window_t));
      % 3) return Root of mean squared signal
      c.data = sqrt(c.data);
    end

    % discard suffixes from contfn/contfilt, we'll provide our own
    c.chanlabels = oldcl; 
    
   otherwise
    error('unsupported envelope ''method''');
  end

  % create new chanlabels
  
  if ~isempty(c.chanlabels) && ~a.nosuffix
    for k = 1:nchans,
      c.chanlabels{k} = [c.chanlabels{k} suffix];
    end
  end
  
  % update data range
  c = contdatarange(c);
  
  % data integrity check
  contcheck(c);