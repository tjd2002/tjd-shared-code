function c = imcont(varargin)
% IMCONT make contdata struct from timeseries/mwl eeg struct/.eeg file
%
% possible inputs and their implied options:
%  
%  'eegfile'/'mwlIOeegfh' = name/filehandle (from mwlopen) of .eeg file
%  -'chans' = vector of channels to use (required)
%  -'timeunits' = 'timestamps' as returned from mwleegfile/load
%  -'data' = data (m x n channels)
%    -if either data or timestamp not provided, load both from file
%    -'timeunits' defaults to 'timestamps'
%
%  'eeg' = struct from eeg2mat
%   -'timeunits' defaults to 'seconds'
%   -'dataunits' is derived from eeg.info, converted to mV
%   -can't select chans, just use all
%  
%  'neuralynxCSC' = struct from NlxLoadCSC
%    (must have 'samptimes', 'samples', 'timeunits', 'dataunits'  fields)
%   -'timeunits' specified by struct
%   -'dataunits' specified by struct
%   -'invert' derived from ncs.info.header.InputInverted 
%   -'chanlabels' specified from ncs.info.header.AcqEntName
%
%  'timestamp' = times in seconds
%  'timestamp_ends' = time of first/last sample
%  'data' = numeric data (m x n channels)
%    -samplerate calcd from timestamps
%    -'timeunits' defaults to 'seconds'
%
%  Common options:
  
%  'invert' : Whether to invert the sign of the data. (default: true
%             whenever 'eegfile', 'mwlIOeegfh', or 'eeg' struct is
%             provided. false if only 'data'/'timestamp' provided.)
%  'resample': resamples data (using 'resample' or 'decimate' as
%              appropriate) by given factor (use 1/integer to decimate)
%  'timewin': range of times to select (selected range will enclose request)
%  'name': a name for the contstruct (default based on filename, chans)
%  'convertdatafun' : defaults to '@single'; use [] to not convert
%  'chanvals': for 2-d data, values associated with each column (slopes,
%              e.g.) for plotting as images/surfs;
%  'chanlabels': text labels for channels
%  'dataunits': defaults to '' (attempts to convert to mV if possible)
%  'timeunits': 'timestamps', 'microseconds' or 'seconds'; default depends on inputs
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
% $$$      -datarange (advisory limits of data -- useful for plotting, but not guaranteed to be real limits of data)
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

  mwl_ts_per_sec = 1e4; % 10000 timestamps/second
  
  a = struct(...
      'eegfile','',...
      'mwlIOeegfh', {{}},... % filehandle from mwlopen;
      'epos', [],... % e.pos struct
      'eposflds', {{}},... % field names in an e.pos struct
      'chans', [],...
      'invert', [],...
      'resample', [],...
      'timewin', [],...
      'eeg', {{}},...
      'neuralynxCSC', [],...
      'ncsfiles', {{}},...
      'name', '',...
      'data', [],...
      'timestamp', [],...
      'timestamp_ends', [],...
      'timeunits', '',...
      'convertdatafun', @single,... 
      'chanvals', [],...
      'chanlabels', {{}},...
      'dataunits', '' );
  
  a = parseArgsLite(varargin,a);

  if isempty(a.invert)
    invert = (~isempty(a.mwlIOeegfh) ||...
              ~isempty(a.eegfile) ||...
              ~isempty(a.eeg));
  else
    invert = a.invert;
  end
    
  c = struct('name', [],...
             'data', [],...
             'chanvals', [],...
             'chanlabels', [],...
             'samplerate',[],...
             'tstart',[], 'tend',[],...
             'datarange',[],...
             'units', '',...
             'nbad_start',0,...
             'nbad_end',0);

  % pass through name, if given
  c.name = a.name;
  
  % string for default names of structs
  chanstr = '';
  % chans[1:4] should look like '1234'
  for k = 1:size(a.chans,2)
    chanstr = [chanstr int2str(a.chans(k))];
  end
  
      
  %%%%%%%%%%%%%%%%%%%%%
  % 'eeg' struct from mwlIO
  
  % use eegfile, if provided
  eegfh = [];
  if ~isempty(a.eegfile),
    eegfh = mwlopen(a.eegfile);
  elseif ~isempty(a.mwlIOeegfh),
    eegfh = a.mwlIOeegfh;
  end
  
  if ~isempty(eegfh);
    
    % come up with a default name
    if isempty(a.name)
      c.name = [eegfh.filename '_ch' chanstr];
    end
    
    timeunits = 'timestamps';

    % # of datapts/timestamp
    recsize = eegfh.nsamples;
    nrecords = eegfh.nrecords;
    nchannels = eegfh.nchannels;
        
    if ~isempty(a.data) && ~isempty(a.timestamp);
      if size(a.data,1) ~= (recsize * size(a.timestamp,1)),
        error('data size inconsistent with timestamps/recordsize');
      elseif size(a.data,2) ~= length(a.chans),
        error('data must have same # of columns as ''chans'' entries');
      end
      
      c.data = a.data;
      timestamp = a.timestamp;
    
    else % we have to load the data;
      
      disp('loading all channels...');
      
      if ~isempty(a.timewin),
        % we are going to prune by timewin later, but this
        % saves us from loading too much excess data.
        
        % get actual timestamps, since there can be gaps
        buff_ts = load(eegfh, 'timestamp');
        buff_ts = buff_ts.timestamp;
        buff_t = double(buff_ts) ./ mwl_ts_per_sec;
        
        loadbuffs = {find(buff_t < a.timewin(1) , 1, 'last'),...
                     find(buff_t > a.timewin(2) , 1, 'first')};
        
        if isempty(loadbuffs{1}),
          loadbuffs{1} = 1;
        end
        
        if isempty(loadbuffs{2}),
          loadbuffs{2} = nrecords;
        end
        
        loadbuffs = [loadbuffs{:}];
        
        if diff(loadbuffs) == 1;
          % need at least 2 timestamps to calc samplerate, etc...
          loadbuffs(2) = loadbuffs(1) + 1;
        end

        t_per_buff = recsize ./ samplerate_reported;
        if max(diff(buff_t(loadbuffs(1):loadbuffs(2)))) > ...
              (1.5 * t_per_buff),
          warning('requested eeg range contains gaps (recording stopped?)'); %#ok
        end

        % load is zero-indexed 
        loadbuffs = loadbuffs-1;
        
        cont_tmp = load(eegfh,'all', loadbuffs(1):loadbuffs(2));
      else
        cont_tmp = load(eegfh,'all', []);
      end
      
      disp('Done...');
      
      if ~isempty(a.convertdatafun) && ~strcmp(func2str(a.convertdatafun), class(a.data)),
          disp(['converting to ' func2str(a.convertdatafun) '...']); % usu. single
          % select appropriate channels, put in argsstruct
          c.data = a.convertdatafun(cont_tmp.data(a.chans,:)');
          disp('Done...');
      end
      
      timestamp = cont_tmp.timestamp(:);
      clear cont_tmp;
      
    end
    
    if ~isempty(a.dataunits)
      c.units = a.dataunits;
    else
      
      % it's in ADCunits; convert to mV;
      gains = getHeaderGains(eegfh);
      gains = gains(a.chans);
      
      % creates a vector of conv factors according to gains
      adunits_to_mv_f = ...
          1/4095 .* ... % ADCrange/ADCunits (-2048 -> +2047)
          20 .* ... % ADCvolts/ADCrange (-10 -> +10)
          1./gains .*... % volts/ADCvolts (vector)
          1000; % mv/volt
      
      % when gain is 0, conversion factor should be 0, not inf
      adunits_to_mv_f(isinf(adunits_to_mv_f)) = 0;
      
      for k = 1:length(a.chans),
        c.data(:,k) = c.data(:,k) .* adunits_to_mv_f(k);
      end
      
      % use full range of AD card as data range
      c.datarange = [-2048 .* adunits_to_mv_f(k),...
                     -2047 .* adunits_to_mv_f(k)];
      c.units = 'mv';
    end
    
    
  %%%%%%%%%%%%%%%%%%%%%
  % struct from eeg2mat
  
  elseif ~isempty(a.eeg), % struct from 'eeg2mat'
    
    % come up with a default name
    if isempty(a.name)
      c.name = [a.eeg.info.filename '_ch' chanstr];
    end

    % eeg struct as input
    timestamp = a.eeg.timestamp;
    timeunits = 'seconds';

    % # of datapts/timestamp
    recsize = 1; 
    
    if strcmp(a.eeg.info.units, 'ADC');
      %convert to mV
      c.units = 'mV';
      
      % a vector of conversion factors per channel
      adunits_to_mv_f = ...
          1/diff(a.eeg.info.ADCrange) .* ... % ADCrange/ADCUnits (-2048 -> +2047)
          diff(a.eeg.info.inputrange) .* ... % ADCvolts/ADCrange (-10 -> +10)
          1/a.eeg.info.gain .* ... % volts/ADCvolts (vector)
          1000; % mv/volt
          
      for k = 1:length(a.eeg.info.gain);
        c.data(:,k) = a.eeg.data(:,k) .* adunits_to_mv_f(k);
      end
      
      c.datarange = [a.eeg.info.ADCrange(1) .* adunits_to_mv_f(k),...
                     a.eeg.info.ADCrange(2) .* adunits_to_mv_f(k)];

      % old code for this line:
% $$$       % good enough: if inputrange doesn't span 0, and there are
% $$$       % different gains, this could be wrong, but who cares?
% $$$       c.datarange = a.eeg.info.ADCrange .* max(adunits_to_mv_f);
      

      
    else
      c.data = a.eeg.data;
      c.datarange = a.eeg.info.ADCrange;
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    % e.pos struct
  elseif ~isempty(a.epos),
    
    recsize = 1;
    
    if isempty(a.eposflds),
      error('if ''epos'' provided, must specify ''eposflds'' to import');
    end
    
    if ischar(a.eposflds),
      a.eposflds = {a.eposflds};
    end
    
    c.data = getepd(a.epos, a.eposflds{:});
    timestamp = getepd(a.epos, 'time');
    
    c.chanlabels = a.eposflds;
    
    if isempty(a.timeunits),
      timeunits = 'seconds';
    else
      timeunits = a.timeunits;
    end

  %%%%%%%%%%%%%%%%%%%%%
  % neuralynx CSC data as imported by Tom's NlxLoadCSC.m
  elseif ~isempty(a.neuralynxCSC)
      
      if ~isempty(a.chanlabels),
          cl = a.chanlabels;
      else
          cl = a.neuralynxCSC.info.header.AcqEntName;
      end
      
      %% Error if user tries to pass in data already provided by the struct
      if ~isempty(a.dataunits),
          error('Data units already specified in neuralynxCSC struct.');
      end

      if ~isempty(a.data),
          error('Data already specified in neuralynxCSC struct.');
      end
      
      if ~isempty(a.timestamp) || ~isempty(a.timestamp_ends) || ~isempty(a.timeunits),
          error('Time (with units) already specified in neuralynxCSC struct.');
      end

      if ~isempty(a.invert)
          error('Data invert status already specified in neuralynxCSC struct.');
      end

      
      %% All other parameters are passed through to the next imcont call
      
      % call imcont recursively with raw timestamps/data
      c = imcont(...
          'timestamp', a.neuralynxCSC.samptimes, ...
          'data', a.neuralynxCSC.samples, ...
          'invert', a.neuralynxCSC.info.header.InputInverted, ...
          'chanlabels', cl,...
          ...
          'chans', a.chans,...
          'resample', a.resample,...
          'timewin', a.timewin,...
          'name', a.name,...
          'chanvals', a.chanvals);
      
      % save Neuralynx-specific info
      c.nlx_info = a.neuralynxCSC.info;

      return;
    
  %%%%%%%%%%%%%%%%%%%%%
  % raw data/timestamp provided by user
  
  else % non-struct inputs (raw 'timestamp'/'data')
    
    recsize = 1;
    
    if isempty(a.data),
      error('must provide ''data''');
    end
    
    c.data = a.data;
    c.units = a.dataunits;

    
    if ~isempty(a.timestamp),
      timestamp = a.timestamp(:);
      timeunits = 'seconds';
    elseif ~isempty(a.timestamp_ends);
      timestamp = linspace(a.timestamp_ends(1), a.timestamp_ends(2), ...
                           size(c.data,1))';
      timeunits = 'seconds';
    else
      timestamp = (1:size(c.data,1))';
      timeunits = 'timestamps';
      warning('no timestamps provided; using integer series');
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%
  % common to all inputs
  
  if size(c.data,1) ~= (recsize * size(timestamp,1)),
    error('data size inconsistent with timestamps/recordsize');
  elseif ~isempty(a.chans) && (size(c.data,2) ~= length(a.chans)),
    error('data must have same # of columns as ''chans'' entries');
  end
  
  [nsamps ncols] = size(c.data); %#ok
  
  % invert sign of data if required
  if invert,
    c.data = -c.data;
  end
  
  % user override of timestamp unitis 
  if ~isempty(a.timeunits);
    timeunits = a.timeunits;
  end
  
  %%% convert timestamps to seconds, if nec.
  switch timeunits
      case 'seconds'
          % no conversion
          if timestamp(end) > 50 * mwl_ts_per_sec;
              warning(['timestamps reported to be in seconds, may be in timestamp ' ...
                  'units or microseconds']);
          end

      case 'microseconds'
          timestamp = double(timestamp) ./ 1e6;

      case 'timestamps'
          timestamp = double(timestamp) ./ mwl_ts_per_sec;

      otherwise
          error('unrecognized ''timeunits'' string');
  end
  
  % pass through chanvals, if given
  if ~isempty(a.chanvals) && length(a.chanvals) ~= ncols,
    error('''chanvals'' must have as many entries as data columns');
  else
    c.chanvals = a.chanvals(:)';
  end
  
  % pass through chanlabels, if given
  if ~isempty(a.chanlabels),
    if ischar(a.chanlabels),
        a.chanlabels = {a.chanlabels};
    end
    if length(a.chanlabels) ~= ncols,
      error('''chanlabels'' must have as many entries as data columns');
    else
      c.chanlabels = a.chanlabels(:)'; %row
    end
  end
  
  %%% calculate samplerate
  
  % IMPORTANT: computed, not reported samplerate, so that we can
  % synthesize new timestamps that sync up with other timestamp-locked
  % events (like spikes): elapsed sample periods / elapsed time  

  %%% crop data to timewin, if requested
  %%% (if from eeg file, already cropped on load)
  if isempty(eegfh) && ~isempty(a.timewin);
 
    % time of first data pt
    sampstart = find(timestamp<a.timewin(1),1,'last');
    if ~isempty(sampstart)
      samprange(1) = sampstart;
    else
      samprange(1) = 1;
    end

    sampend = find(timestamp>a.timewin(2),1,'first');
    if ~isempty(sampend)
      samprange(2) = sampend;
    else
      samprange(2) = size(timestamp,1);
    end
        
    samprange = max(1, samprange);
    samprange = min(size(c.data,1), samprange);
    
    c.tstart = timestamp(samprange(1));
    c.tend = timestamp(samprange(2));
    
    % # of sample periods elapsed / time elapsed
    c.samplerate = diff(samprange)*recsize / (c.tend-c.tstart);
    
    % crop data/timestamp
    c.data = c.data(samprange(1):samprange(2),:);
    timestamp = timestamp(samprange(1):samprange(2));
    
    % test/warn timewin
    if c.tstart > a.timewin(1) || c.tend < a.timewin(2),
      warning('requested timewin out of data range, used available data');
    end
  
  else
    
    % # of sample periods elapsed / time elapsed
    c.samplerate = (size(timestamp,1)-1)*recsize / diff(timestamp([1 end]));
    
    % time of first data pt
    c.tstart = timestamp(1);
    
    % time of last data pt
    c.tend = timestamp(end) + (recsize-1)/c.samplerate;
  end
  
  
  %%% calculate max error between computed/exact timestamps
  %
  % compare actual timestamps to 'synthetic' even sampling
  %
  % [could be more memory efficient using a cumsum of real ts and a
  % counter for real ts, but !whatever!]
  ts_syn = c.tstart:(recsize/c.samplerate):c.tend;
  c.max_tserr = max(abs(timestamp(:) - ts_syn(:)));
  
  % there's no reason the error shouldn't be 0, but 50% of a sample is a
  % natural cutoff/confidence interval. Note that this error is the max of
  % what the synthesizing of the timestamps **adds** to the timing
  % uncertainty already in the signal.(for eeg, this should be only 0.5
  % of a *timestamp interval* = 50usec), so really no biggie.
  
  if c.max_tserr > 0.5*(1/c.samplerate),
    warning('imcont:ts_error', ['max error in computed timestamp will be greater than 50% of ' ...
             'samplerate resolution!!!--try converting as subsets of data?']);
  end
  
  %%% resample, if requested
  
  if ~isempty(a.resample),
    disp('resampling...');
    c = contresamp(c,'resample',a.resample);
    disp('Done...');
  end
    
    
  %%% useful when plotting; nice to precompute, even though it's slow
  c = contdatarange(c);
  
  %%% convert datatype last-ish
  if ~isempty(a.convertdatafun) && ~strcmp(func2str(a.convertdatafun), class(a.data)),
      disp(['converting to ' func2str(a.convertdatafun) '...']); % usu. single
    c.data = a.convertdatafun(c.data);
  end
  
