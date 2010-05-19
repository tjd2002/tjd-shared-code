function c = imconthist(varargin)
% IMCONTHIST make contdata struct from event list (histogram)
%
% pass in evtimes, values, can threshold on values
%
% possible inputs:
%  
%  'cls': e.cl structure array. If provided, we will extract evtimes, 
%         chanlabels, and optionally data (see next option)
%  'cls_thresh_fld': optional named field to use for data for thresholding
%  'cls_merge': if true, merge all events in all clusters into a single histogram
%
%  'evtimes' = m x 1 numeric array of event times.
%  'timeunits' defaults to 'seconds'
%  'data' (optional) (m x 1 numeric data to use for thresholding)
%
%  'histbinsize' size of hist bins, in timeunits
%
%  if 'data' or 'cls_thresh_fld' is provided, we can select which events
%  to include in the histogram.
%     'thresh_val', only 
%     'thresh_fn' defautls to '@gt', greater-than, >
%
%  'timewin': range of evtimes to select (selected range will enclose request)
%  'name': a name for the contstruct (default based on filename, chans)
%  'convertdatafun' : defaults to '@single'; use [] to not convert
%  'chanlabels': for 2-d data, values associated with each column (slopes,
%              e.g.) for plotting as images/surfs;
%  'convert_freq': if true, report hist as frequency in Hz, rather than event
%          counts/bin. (default false)
%
% output struct:
%
% $$$    -contdata struct:
% $$$      -name
% $$$      -data (anything numeric)
% $$$      -chanvals (for 2-D data, values associated with each column)
% $$$      -chanlabels (for 2-D data, labels for each column) cell array of strings
% $$$      -samplerate (derived, not reported; in samples/second)
% $$$      -tstart,tend (time of first/last sample; in seconds)
% $$$      -datarange (advisory limits of data -- useful for plotting, but
% $$$       not guaranteed to be real limits of data)
% $$$      -units (descriptive string)
% $$$      -nbad_start and nbad_end: samples that are unreliable due to filtering
% $$$
% $$$      -[sourcefile/channel]
% $$$      -[prefilt - analog filtering by amplifier]
% $$$      -[filters{m} - actual filters applied (using filtfilt)]
%
% todo:
% -input checking before all the loading/parsing/resampling, to save user
% hassle, esp chanvals and chanlabels

  ts_per_sec = 1e4; % 10000 timestamps/second
  
  a = struct(...
      'cls', [],...
      'cls_thresh_fld', '',...
      'cls_merge', false,...
      'evtimes', [],...
      'data', [],...
      'histbinsize', [],...
      'thresh_val',[],...
      'thresh_fn',@gt,...
      'timewin', [],...
      'name', '',...
      'timeunits', 'seconds',...
      'convertdatafun', @single,... 
      'chanlabels', {{}},...
      'convert_freq', false);
  
  a = parseArgsLite(varargin,a);
  
  c = struct('name', [],...
             'data', [],...
             'chanvals', [],...
             'chanlabels', [],...
             'samplerate',[],...
             'tstart',[], 'tend',[],...
             'datarange',[],...
             'units', '',...
             'nbad_start',0,...
             'nbad_end',0,...
             'max_tserr', []);

  % require histbinsize
  if isempty(a.histbinsize),
    error('must provide ''histbinsize''');
  end
  
  % pass through name, if given
  c.name = a.name;

  if length(a.cls) < 2,
    a.cls_merge = false;
  end
  
  % parse timeunits
  switch(a.timeunits)
   case 'timestamps', 
    time_in_ts = true;
   case 'seconds', 
    time_in_ts = false;
   otherwise
    error('unrecognized ''timeunits'' string');
  end

  % pass through chanlabels, if given
  if ~isempty(a.chanlabels) && length(a.chanlabels) ~= ncols,
    error('''chanlabels'' must have as many entries as data columns');
  else
    c.chanlabels = a.chanlabels(:)'; %row
  end
  
  evtimes = {};
  datas = {};
  
  if ~isempty(a.cls),
    for k = 1:length(a.cls)
      evtimes{k} = getepd(a.cls(k),'time');
      if ~isempty(a.cls_thresh_fld),
        datas{k} = getepd(a.cls(k),a.cls_thresh_fld);
      end
    
      if isempty(a.chanlabels),
        if a.cls_merge
          % only one channel, label it
          label = strvcat(a.cls.name); %#ok
          label = [label repmat('_',size(label,1),1)]';
          % leave off last char b/c it's a '_'
          c.chanlabels{1} = ['merge_' label(1:end-1)];
        else
          c.chanlabels{k} =  a.cls(k).name;
        end
      end
    end
    
    if a.cls_merge,
      evtimes = {vertcat(evtimes{:})};
      datas = {vertcat(datas{:})};
    end
    
    
  else
    % no cls, get evtimes, data from args
    
    if ~iscell(a.evtimes),
      evtimes = {a.evtimes};
    else
      evtimes = a.evtimes;
    end
    
    if ~iscell(a.data),
      datas = {a.data};
    end
  end

  
  % if no timewin provided, calculate
  if isempty(a.timewin) || any(isinf(a.timewin))
    if isempty(evtimes),
      error('no/inf timewin and no events; can''t make histogram');
    else
      warning('no timewin, using first/last event as window');
      alltimes = (vertcat(evtimes{:}));
      a.timewin = [ min(alltimes) max(alltimes) ];
    end
  end
  
  
  % make histbins
  histbin_edges = a.timewin(1):a.histbinsize:a.timewin(2);
  
  % loop over cell array of evtimes
  for k = 1:length(evtimes)
    evt = evtimes{k};

    % select above-threshold events
    if ~isempty(a.thresh_val),
      evti = a.thresh_fn(datas{k}, a.thresh_val);
      evt = evt(evti);
    end    

    % get event histogram (-Inf is hack to deal with empty evt)
    evt_hist = histc([-Inf; evt(:)], histbin_edges);

    % discard last bin (values equal to highest edge, zero-width bin)
    evt_hist = evt_hist(1:end-1);
    
    % get average firing rate in Hz by scaling by: # of cells contributing;
    % sampling rate
    if a.convert_freq,
      c.units = 'Hz';
      evt_hist = evt_hist./a.histbinsize;

      % if timeunits are timestamps, extra conversion to get Hz
      if time_in_ts
        spks_hist = spks_hist .* ts_per_sec;
      end
    else
      c.units = 'counts/bin';
    end
    
    % if merging multiple clusters, get average f.r.per cluster
    if a.cls_merge,
      evt_hist = evt_hist./length(a.cls);
    end

    % store each 
    c.data(:,k) = evt_hist;
  end
  
  % bin centers are appropriate place to plot each value
  c.tstart = histbin_edges(1)+ (0.5*a.histbinsize);
  c.tend = histbin_edges(end)- (0.5*a.histbinsize);

  % convert from timestamps if requested
  if time_in_ts,
    c.tstart = double(c.tstart) ./ ts_per_sec;
    c.tend = double(c.tend) ./ ts_per_sec;
  end

  % samplerate is exact, no ts error introduced by histc
  c.samplerate = 1/a.histbinsize;
    
  %%% useful when plotting; nice to precompute, even though it's slow
  c = contdatarange(c);
  
  % convert datatype last-ish
  if ~isempty(a.convertdatafun),
    c.data = a.convertdatafun(c.data);
  end
  
  
  if ~time_in_ts && (c.tend - c.tstart) > (50 * ts_per_sec);
    warning(['timestamps reported to be in seconds, may be in timestamp ' ...
             'units']);
  end
  
