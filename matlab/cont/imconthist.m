function c = imconthist(varargin)
% IMCONTHIST - create rate histogram cont structs from lists of event times 
%
%  c = imconthist([name/value pairs])
%
% By default, creates a cont struct with one channel per event list
% provided, with rate in Hz.
%
% Inputs: (* means required, -> indicates default value)
%
%  Basic Inputs:
%  *'evtimes' - numeric vector of event times (or cell array of vectors)
%  *'timewin': range of time over which to calculate histogram, in 'timeunits'.
%       (If you don't provide a timewin, imconthist will treat the first and
%       last event times as the timewin; this could introduce bias into
%       edge bins.)
%  *'histbinsize' - size of hist bins, in 'timeunits'.
%   'timeunits' - (->'seconds', 'microseconds', 'MWL_timestamps',
%       'NLX_timestamps')
%   'convert_freq' - if true, return rate as events/s, rather than event
%       counts/bin (default true).
%   'merge' - if true, merge all events provided into a single histogram
%       (default false).
%   'mergenorm' - if true, normalize (divide) histogram by number of input
%       lists, giving average rate; if false, sum events in each time bin
%       (default true).
%
%  Inputs controlling filtering of events based data associated with each
%  event:
%   'data' - data associated with each spike (vector or cell array of
%       vectors same size as 'evtimes')
%   'thresh_val' - threshold to apply to each data
%   'thresh_fn' - thresholding function (defaults to '@gt', greater-than)
%
%   'name' - a name for the contstruct (default based on filename)
%   'chanlabels' - cell array of strings to label each new channel
%   'convertdatafun' - defaults to '@single'; use [] to not convert
%   'suffix' - string to append to computed chanlabels (default 'rate')
%
%  Specialized inputs:
%   Instead of evtimes/data, can provide an array of cluster (.cl)
%   structures containing named data fields:
%    'cls': e.cl structure array. If provided, we will extract evtimes, 
%         chanlabels, and optionally data (see next option)
%    'cls_thresh_fld': optional named field to use for data for thresholding
%
%
% Outputs:
%  c: cont struct ('help imcont' for details)
%   -data is in events/timeunit (so, by default, Hz)
%
%
% Example: From a list of spike times for 3 electrodes, (in variables t1, t2, and
%     t3, seconds), and corresponding lists of peak amplitude for each spike (in
%     pks1-pks3, mV), calculate average firing rate per electrode in Hz, in 25ms
%     bins, but only including spikes that reach at least 100uV, and only
%     for the period 1000-2000s.
%
%   cont_mua = imconthist('evtimes', {t1 t2 t3},...
%                         'timewin', [1000 2000],...
%                         'histbinsize', 0.025,...
%                         'convert_freq', true,...
%                         'data', {pks1 pks2 pks3},...
%                         'thresh_val', 0.100,...
%                         'merge', true,...
%                         'mergenorm', true,...
%                         'chanlabels', {'MUA rate, cells 1-3'}) 
%
% Tom Davidson <tjd@stanford.edu> 2003-2010
  
% TODO:
%  -change histbinsize, timewin, to be in seconds, regardless of
%  'timeunits'. (change name of timeunits to 'evtimeunits' to reflect
%  change).

%% constants
  AD_ts_per_sec = 1e4; % MWL's 'AD' timestamps
  NLX_ts_per_sec = 1e6; % Neuralynx timestamps
  
  %% parse inputs
  a = struct(...
      'evtimes', [],...
      'data', [],...
      'cls', [],...
      'cls_thresh_fld', '',...
      'merge', false,...
      'mergenorm', true,...
      'histbinsize', [],...
      'thresh_val',[],...
      'thresh_fn',@gt,...
      'timewin', [],...
      'name', '',...
      'timeunits', 'seconds',...
      'convertdatafun', @single,... 
      'chanlabels', {{}},...
      'suffix', 'rate',...
      'convert_freq', true);
  
  a = parseArgsLite(varargin,a);

  c = mkcdat('name', a.name);
  
  % require histbinsize
  if isempty(a.histbinsize),
    error('must provide ''histbinsize''');
  end
  
  % only one input method at a time, please
  if ~isempty(a.cls) && ~isempty(a.evtimes)
    error('Can only provide one of ''cls'', and ''evtimes''');
  end
  
  % parse timeunits
  switch lower(a.timeunits)
   case 'MWL_timestamps', 
    ts_per_sec = AD_ts_per_sec;
   case {'NLX_timestamps'};
    ts_per_sec = NLX_ts_per_sec;
   case {'microseconds'}
    ts_per_sec = 1e6;
   case 'seconds', 
    ts_per_sec = 1;
   otherwise
    error('unrecognized ''timeunits'' string');
  end

  evtimes = {};
  datas = {};
  
  %% get evtimes and datas from cls
  if ~isempty(a.cls),
    for k = 1:length(a.cls)
      evtimes{k} = getepd(a.cls(k),'time');
      if ~isempty(a.cls_thresh_fld),
        datas{k} = getepd(a.cls(k),a.cls_thresh_fld);
        usedata = true;
      else
        usedata = false;
      end
    end
    
  else

    %% no cls, get evtimes, data from args
    
    % convert to cell arrays if nec
    if ~iscell(a.evtimes),
      evtimes = {a.evtimes};
    else
      evtimes = a.evtimes;
    end
    
    % are we going to select events using data?
    usedata = ~isempty(a.data);   
    if usedata
      if ~iscell(a.data),
        datas = {a.data};
      else
        datas = a.data;
      end

      % make sure we have matching numbers of time/data inputs
      if numel(evtimes) ~= numel(datas),
        error(['If ''data'' provided, must provide one vector for each list     ' ...
               'of ''evtimes''.']);
      end
      
    end
  end
  
  % make sure evtimes/data are equal-sized column vectors
  for k = 1:numel(evtimes)

    if ~isvector(evtimes{k})
      error('Times must be provided as vectors');
    else
      % convert to column vector
      evtimes{k} = evtimes{k}(:);
    end
    
    if usedata
      if any(size(datas{k} ~= size(evtimes{k})))
        error(['''data'' lists must be same size as ''evtimes'' lists ' ...
               '(#%d)']);
      end
      if ~isvector(datas{k})
        error('Times must be provided as vectors');
      else
        datas{k} = datas{k}(:);
      end
    end
  end

  % how many inputs (before merging)?
  ninputs = numel(evtimes);
  
  % merge times / datas across lists if requested
  if a.merge,
    evtimes = {vertcat(evtimes{:})};
    datas = {vertcat(datas{:})};
  end
  
  % how many histograms?
  nhists = numel(evtimes);
  
  % if no timewin provided, calculate
  if isempty(a.timewin) || any(isinf(a.timewin))

    % find min/max event time across all inputs (memory efficient)
    t_min = Inf; t_max = -Inf;
    for k = 1:nhists
      t_min = min(t_min, min(evtimes{k}));
      t_max = max(t_max, max(evtimes{k}));
    end
    
    if isinf(t_min),
      % there were no events in any list
      error('No ''timewin'' and no events; can''t make histogram');
    end

    warning('No ''timewin'' provided, so using first/last event as window');
    timewin = [t_min t_max];
  
  else
    % use provided timewin
    timewin = a.timewin;
  end

  if diff(timewin)<0,
    error('Timewin must be non-decreasing: Timewin = [%d %d]', timewin(1), ...
          timewin(2));
  end
  if diff(timewin)<a.histbinsize,
    error('Timewin must be at least 1 hist bin wide: Timewin = [%d %d]', timewin(1), ...
          timewin(2));
  end
  
  % make histbins
  histbin_edges = timewin(1):a.histbinsize:timewin(2);
  
  % loop over cell array of evtimes
  for k = 1:nhists
    evt = evtimes{k};

    % select events meeting threshold
    if ~isempty(a.thresh_val),
      evti = a.thresh_fn(datas{k}, a.thresh_val);
      evt = evt(evti);
    end    

    % get event histogram (Inf is hack to deal with empty evt list)
    evt_hist = histc([evt(:); Inf], histbin_edges);

    % discard last bin (values equal to highest edge, zero-width bin)
    evt_hist = evt_hist(1:end-1);
    
    % get average firing rate in Hz by scaling by: # of cells contributing;
    % sampling rate
    if a.convert_freq,
      c.units = 'Hz';
      evt_hist = evt_hist./a.histbinsize;

      % if timeunits are timestamps, extra conversion to get Hz
      if ts_per_sec~=1
        spks_hist = spks_hist .* ts_per_sec;
      end
    else
      c.units = 'counts/bin';
    end
    
    % if merging multiple clusters, get average f.r.per cluster
    if a.merge && a.mergenorm
      evt_hist = evt_hist./ninputs;
    end

    % store each 
    c.data(:,k) = evt_hist;
  end
  
  % bin centers are appropriate place to plot each value
  c.tstart = histbin_edges(1)+ (0.5*a.histbinsize);
  c.tend = histbin_edges(end)- (0.5*a.histbinsize);

  % convert from timestamps if requested
  if ts_per_sec ~= 1,
    c.tstart = double(c.tstart) ./ ts_per_sec;
    c.tend = double(c.tend) ./ ts_per_sec;
  end

  % samplerate is exact, no ts error introduced by histc
  c.samplerate = 1/a.histbinsize;

  %% Make some nice chanlabels
      
  % pass through chanlabels, if given
  if ~isempty(a.chanlabels)
    if numel(a.chanlabels) ~= nhists,
      error('''chanlabels'' must have as many entries as evtimes/cls');
    else
      c.chanlabels = a.chanlabels(:)'; % make a row
    end
  else
    % no labels provided, can we make them from the .cls?
    if ~isempty(a.cls)
      if a.merge
        % only one channel, label it
        label = strvcat(a.cls.name); %#ok
        label = [label repmat('_',size(label,1),1)]';
        % leave off last char b/c it's a '_'
        c.chanlabels{1} = ['merged_' label(1:end-1)];
      else
        c.chanlabels{k} =  a.cls(k).name;
      end
    else
      % no chanlabels, no .cls, just number the channels
      if a.merge
        c.chanlabels{1} = 'merged';
      else 
        for k = 1:nhists
          c.chanlabels{k} = num2str(k);
        end
      end
    end
    for k = 1:nhists
      c.chanlabels{k} = [c.chanlabels{k} '_' a.suffix];
    end
  end

  % update data range
  c = contdatarange(c);
  
  % convert datatype last-ish
  if ~isempty(a.convertdatafun),
    c.data = a.convertdatafun(c.data);
  end
  
  if ts_per_sec==1 && (c.tend - c.tstart) > (50 * ts_per_sec);
    warning(['timestamps reported to be in seconds, but may be in timestamp ' ...
             'units (range: %d - %d'], c.tstart, c.tend);
  end
  
  % data integrity check
  contcheck(c);
