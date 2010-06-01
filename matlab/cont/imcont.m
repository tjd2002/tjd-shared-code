function [c timestamp ts_syn] = imcont(varargin)
% IMCONT - create cont data struct from timeseries data
%
% Creates a 'cont' data struct from timeseries data (fields documented in
% 'outputs', below). Multiple channels with the same time base can be stored
% in the same cont struct. Cont structs can be operated on by the cont*
% family of functions (contfilt, contcombine, contresamp, contwin ...), most
% of which produce another cont struct as their output.
%
% NB there are other functions to import data into cont structs:
%  imconthist - create a rate histogram from a list of event times.
%  imcontsegs - create a cont struct from discontinuous chunks of
%      timeseries data.
%
% Inputs: (* = required, -> indicates a default value)
%  You can provide timeseries data in a number of formats, including
%  data in variables, or paths to files. Where available, sampling rate,
%  data units, and other signal metadata (e.g. channel name) is extracted
%  from the data files.
%
%   Raw Data:
%   *'data' - numeric data (m x n channels)
%    'timestamp' - time in seconds at each data point
%    'timestamp_ends' - time of first/last sample
%    Effect on other inputs:
%     'timeunits' - defaults to 'seconds'
%     'invert' defaults to false
%
%   Neuralynx .csc data struct, as imported by Tom's NlxLoadCSC function
%   *'neuralynxCSC' - struct from NlxLoadCSC
%    Effect on other inputs:
%     'timeunits' = specified by struct
%     'dataunits' = specified by struct
%     'invert' = specified by header's 'InputInverted' field
%     'chanlabels' - defaults to header's 'AcqEntName' field
%
%   Wilson Lab (MIT) .eeg file - Note: depends on Fabian's mwlIO library
%   *'eegfile'/'mwlIOeegfh' = filename/filehandle (from mwlopen) of .eeg file
%    'chans' - vector of channels to use (import all by default)
%    Effect on other inputs:
%     'timeunits' = 'MWL_timestamps'
%     'invert' - (default: true)
%  
%   E.pos struct (as used in Tom's expt struct)
%   *'epos' - the epos struct
%   *'eposflds' - names of epos fields to import as channels (string / cell
%       array of strings)
%    Effect on other inputs:
%     'timeunits' - defaults to 'seconds'
%     'invert' - defaults to false
%
%  Shared Inputs: (may be implied/overwritten depending on import type, above)
%   'invert' - Whether to invert the sign of the data. (see above for defaults)
%   'timewin'- range of times to select (seconds).
%   'name'- a name for the contstruct (default based on filename, chans).
%   'convertdatafun' - defaults to '@single'; use [] for no conversion.
%   'chanvals'- for 2-d data, scalar values associated with each channel.
%   'chanlabels'- text labels for channels.
%   'dataunits'- defaults to '' (attempts to convert to mV if possible).
%   'timeunits'- 'MWL_timestamps', 'microseconds' or 'seconds'.

%   'ts_syn_linmode' - method for fitting a linear samplerate to observed
%       sample times. 'ends' (the default) simply linearizes from first to last
%       observed timestamp. 'regress' performs a least-squares linear
%       regression and should give lower mean error with less bias (i.e. possibly
%       higher max_tserr, but lower mean_tserr)
%
%   'allowable_tserr_samps' - maximum allowable timestamp error, in # of
%       samples (default 0.5).
%
% Output:
%   c - cont struct with fields:
%     'name' - string, cont functions often modify this.
%     'data' - 2-D array of numeric data; 1 row per sample, 1 column per
%         channel. All numeric datatypes supported, defaults to
%         single-precision floating point.
%     'chanlabels' - for multi-channel data, string labels for each column
%     'chanvals' - for 2-D data, scalar values associated with each
%         column (e.g. frequency at each bin for a spectrogram).
%     'samplerate' - derived, not claimed (samples/second, double-precision)
%     'tstart','tend' - time of first/last sample; in seconds.
%     'datarange' - advisory limits of data -- useful for plotting, but
%         not guaranteed to be real limits of data; in data units.
%     'units' - descriptive string giving units of data.
%     'nbad_start' and 'nbad_end' - samples at start/end of data that are
%         unreliable due to filtering edge effects)
%     'max_tserr' - maximum timestamp error introduced across entire signal
%         by the conversion from 1 timestamp/sample, to a fixed
%         samplerate.
%     'mean_tserr' - mean timestamp error, or mean temporal bias.
%
% Example: Import a Neuralynx .csc file, downsample to exactly 1kHz
%
%   cd data/lfpdir
%   csc = NlxLoadCSC('LFP1.ncs')
%   cdat = imcont('neuralynxCSC', csc);
%   cdat = continterp(cdat, 'samplerate', 1000);
%
% If you encounter out of memory errors, try loading in subsets of
% channels, then merging the cont structs with contcombine.
  
% Tom Davidson <tjd@stanford.edu> 2003-2010
  
% TODO:
% -input checking before all the loading/parsing/resampling, to save user
% debugging time, esp chanvals and chanlabels
% -handle calls to NlxLoadCSC, allow 'ncsfiles' argument

  mwl_ts_per_sec = 1e4; % 10,000 timestamps/second
  
  a = struct(...
      'eegfile','',...
      'mwlIOeegfh', {{}},... % filehandle from mwlopen;
      'epos', [],... % e.pos struct
      'eposflds', {{}},... % field names in an e.pos struct
      'chans', [],...
      'invert', [],...
      'resample', [],...
      'timewin', [],...
      'neuralynxCSC', [],...
      'name', '',...
      'data', [],...
      'timestamp', [],...
      'timestamp_ends', [],...
      'timeunits', '',...
      'convertdatafun', @single,... 
      'ts_syn_linmode', 'ends',...
      'allowable_tserr_samps', 0.5,...
      'chanvals', [],...
      'chanlabels', {{}},...
      'dataunits', '' );
  
  a = parseArgsLite(varargin,a);

  if isempty(a.invert)
    % by default, invert is true for MWL files
    invert = (~isempty(a.mwlIOeegfh) ||...
              ~isempty(a.eegfile));
  else
    invert = a.invert;
  end
    
  % call cont data constructor function
  c = mkcdat('name', a.name);
  
  % string for default names of structs
  chanstr = '';
  % chans[1:4] should look like '1234'
  for k = 1:size(a.chans,2)
    chanstr = [chanstr int2str(a.chans(k))];
  end
  

  % Figure out what kind of inputs we're working with
  mode = {};    

  eegfh = [];
  if ~isempty(a.eegfile),
    eegfh = mwlopen(a.eegfile);
    mode = [mode {'mwleeg'}];
  elseif ~isempty(a.mwlIOeegfh),
    eegfh = a.mwlIOeegfh;
    mode = [mode {'mwleeg'}];
  end

  if ~isempty(a.epos)
    mode = [mode {'epos'}];
  end
  
  if ~isempty(a.neuralynxCSC),
    mode = [mode {'nlxCSC'}];
  end

  if ~isempty(a.data) || ~isempty(a.timestamp),
    mode = {'raw'};
  end
  
  switch numel(mode)
   case 0
    error('You must provide data or a file to import');
   case 1
    % ok, 1 mode only
   otherwise
    modestr = strcat(mode, ',');
    error('More than one data import mode implied by inputs: %s.', modestr(1:end-1));
  end
  
  mode = mode{1};

  %% Get c.data and timestamp data from any of the import modes, then do
  %% generic processing below
  switch mode

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Import Wilson lab .eeg file import using Fabian's mwlIO library
   case 'mwleeg'
    
    if ~isempty(a.timeunits)
      error(['Can''t provide timeunits for MWL eeg import--times always ' ...
             'in MWL_timestamp units']);
    end
    timeunits = 'MWL_timestamps';

    % # of datapts/timestamp
    recsize = eegfh.nsamples;
    nrecords = eegfh.nrecords;
    nchannels = eegfh.nchannels;

    %% Load Timestamps, select buffers, and check for gaps
    disp('Loading timestamps...');
    % load actual timestamps, since there can be gaps
    buff_ts = load(eegfh, 'timestamp');
    buff_ts = buff_ts.timestamp;
    buff_t = double(buff_ts) ./ mwl_ts_per_sec;
    
    if ~isempty(a.timewin),
      % we are going to prune by timewin later, but this
      % saves us from loading too much excess data.
      loadbuffs = {find(buff_t < a.timewin(1) , 1, 'last'),...
                   find(buff_t > a.timewin(2) , 1, 'first')};

      if isempty(loadbuffs{1}),
        loadbuffs{1} = 1;
      end
      if isempty(loadbuffs{2}),
        loadbuffs{2} = nrecords;
      end
      loadbuffs = [loadbuffs{:}];
    else
      loadbuffs = [1 nrecords];
    end

    if diff(loadbuffs) == 1;
      % need at least 2 timestamps to calc samplerate
      loadbuffs(2) = loadbuffs(1) + 1;
    end

    % note this is samplerate according to the A->D card's clock, not the master
    % computer's clock. OK for our purposes here, though.
    samplerate_reported = str2double(eegfh.header(2).rate)./nchannels;
    t_per_buff = recsize ./ samplerate_reported;
    if max(diff(buff_t(loadbuffs(1):loadbuffs(2)))) > ...
          (1.1 * t_per_buff),
      warning('Requested eeg range contains gaps (recording stopped?)'); %#ok
    end
    
    disp('Done...');

    %% Load in data
    disp('Loading data...')
    loadbuffs = loadbuffs-1; % 'load' is zero-indexed 
    cont_tmp = load(eegfh,'all', loadbuffs(1):loadbuffs(2));
    
    % select requested channels (data in rows)
    if ~isempty(a.chans)
      c.data = cont_tmp.data(a.chans,:)';
    else
      c.data = cont_tmp.data(:,:)';
    end

    timestamp = cont_tmp.timestamp(:);
    clear cont_tmp;

    disp('Done...');
    
    % Change datatype before converting units since mwlopen returns
    % int16s, which underrun when converted to mV 
    if ~isempty(a.convertdatafun) && ~strcmp(func2str(a.convertdatafun), class(c.data)),
      disp(['converting to ' func2str(a.convertdatafun) '...']); % usu. single
      
      % select appropriate channels, put in argsstruct
      c.data = a.convertdatafun(c.data);
      disp('Done...');
    end

    % OK, can convert units now
    
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
      
      % report full range of AD card as data range
      c.datarange = [-2048 .* adunits_to_mv_f(k),...
                      2047 .* adunits_to_mv_f(k)];
      c.units = 'mV';
    end
    
    % come up with a default cont struct name
    if isempty(c.name)
      c.name = [eegfh.filename '_ch' chanstr];
    end

    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Import e.pos struct from Tom's expt struct    
   case 'epos'
    
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

    c.units = a.dataunits;
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Import Neuralynx CSC data struct (as created by Tom's NlxLoadCSC.m)
   case 'nlxCSC'
    
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

    %% All other parameters are passed through as-is to the next imcont call

    % call imcont recursively with raw timestamps/data options
    c = imcont(...
        'timestamp', a.neuralynxCSC.samptimes, ...
        'data', a.neuralynxCSC.samples, ...
        'invert', a.neuralynxCSC.info.header.InputInverted, ...
        'chanlabels', cl,...
        ...
        'chans', a.chans,...
        'timewin', a.timewin,...
        'name', a.name,...
        'chanvals', a.chanvals,...
        'convertdatafun', a.convertdatafun,... 
        'ts_syn_linmode', a.ts_syn_linmode,...
        'allowable_tserr_samps', a.allowable_tserr_samps);

    % save Neuralynx-specific info
    c.nlx_info = a.neuralynxCSC.info;

    % Don't do additional processing--it was done in call out to imcont
    return;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Import raw data/timestamp lists
   case 'raw'
    
    recsize = 1;
    
    if isempty(a.data),
      error('must provide ''data''');
    end
    
    c.data = a.data;
    c.units = a.dataunits;

    % default to seconds
    if ~isempty(a.timeunits)
      timeunits = a.timeunits;
    else
      timeunits = 'seconds';
    end

    if ~isempty(a.timestamp),
      timestamp = a.timestamp(:);
    elseif ~isempty(a.timestamp_ends);
      timestamp = linspace(a.timestamp_ends(1), a.timestamp_ends(2), ...
                           size(c.data,1))';
    else
      error('No time information provided; can''t create cont struct.');
    end
  end
  
  
  %% Common processing for all input modes
  
  % check reported data size consistency
  if size(c.data,1) ~= (recsize * size(timestamp,1)),
    error('data size inconsistent with timestamps/recordsize');
  elseif ~isempty(a.chans) && (size(c.data,2) ~= length(a.chans)),
    error('data must have same # of columns as ''chans'' entries');
  end
  
  [nsamps ncols] = size(c.data); %#ok
  
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

      case 'MWL_timestamps'
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
  
  %  To save memory and simplify filtering operations, cont structs do not
  % store a timestamp for each data point. Instead, we precisely store the
  % time of the first/last sample and assume even sampling. If the data
  % is not perfectly uniformly sampled, then this procedure introduces
  % timing errors. In the data we analyze, these are due to slight drift
  % between the clock used to drive sampling of our signal, and a master
  % 'timestamp' clock (which is also used to synchronize other systems,
  % such as video recording, digital events, sampling across multiple
  % computers, etc). Both MWL's AD and Neuralynx's Cheetah take this
  % approach.

  %  Put another way, we replace (or model) the observed timestamps with a
  %  linear fit. The errors are usually very small, but before discarding the
  %  true timestamps, we calculate the maximum error introduced by our
  %  model, and save this as a 'worst-case' timing error.

  %% crop data to timewin before synthesizing timestamps so we don't try
  %% to fit error unnecessarily

  % (if from MWL eeg file, already mostly cropped on load)
  if strcmp(mode, 'mwleeg') || isempty(a.timewin);
    % select all records
    tsi_range = [1 numel(timestamp)];
  else 
    % get indexes of selected first and last records
    tsi_range{1} = find(timestamp<a.timewin(1),1,'last');
    if isempty(tsi_range{1})
      tsi_range{1} = 1;
    end

    tsi_range{2} = find(timestamp>a.timewin(2),1,'first');
    if isempty(tsi_range{2})
      tsi_range{2} = size(timestamp,1);
    end
    tsi_range = [tsi_range{:}];
  end
  
  % timestamps of first and last selected records to use
  ts_range = timestamp(tsi_range);
  % index into corresponding samples
  samprange = [((tsi_range(1)-1)*recsize)+1 ... % end of prev record + 1
               tsi_range(2)*recsize]; % end of last record

  % crop to selected data & timestamps
  timestamp = timestamp(tsi_range(1):tsi_range(2));
  c.data = c.data(samprange(1):samprange(2),:);
  
  %%% try to catch large gaps in data (without knowing the samplerate)
  dts = diff(timestamp);
  gapi = find (dts > 1.5*median(dts));
  if ~isempty(gapi)
    error('Timestamp data appears to have gaps: after indexes [%s] (times [%s])',...
            num2str(gapi), num2str(timestamp(gapi)));
  end
  
  %%% calculate linear fit 'synthetic timestamps'
  numts = numel(timestamp);
  switch lower(a.ts_syn_linmode)
   case 'regress'
    % find best fit to real timestamps by linear regression
    tsb = regress(timestamp(:), [ones(numts,1) (1:numts)']);

    % synthesize timestamps from slope/offset
    ts_syn = tsb(1) + (1:numts)'*tsb(2);
    
   case 'robust'
    warning('experimental robust regression!');
    tsb = robustfit((1:numts)', timestamp(:), 'fair', 2);

    % synthesize timestamps from slope/offset
    ts_syn = tsb(1) + (1:numts)'*tsb(2);

   case 'ends'

    % linear spacing between first/last timestamp
    ts_syn = linspace(timestamp(1), timestamp(end), numts)';

   otherwise
    error('Unrecognized ''ts_syn_linmode'': try ''regress'' ''ends''.');
  
  end

  % calculate samplerate, start/end times from ts_syn
  
  % # of sample periods elapsed / time elapsed
  c.samplerate = (size(ts_syn,1)-1)*recsize / diff(ts_syn([1 end]));
  
  % calculate stats on timing error introduced by our linearization procedure
  tserr = timestamp(:) - ts_syn(:);
  [c.max_tserr] = max(abs(tserr));
  [c.mean_tserr] = mean(tserr);
  if c.max_tserr > a.allowable_tserr_samps*(1/c.samplerate),
    error(['max error in computed timestamp (%3.3f ms) will be greater than %.1f%% of ' ...
           'sampling resolution (%3.3f ms) !!!--adjust ''allowable_tserr_samps'' or '...
           'try converting as subsets of data.'],...
          c.max_tserr*1000, a.allowable_tserr_samps*100, 1/c.samplerate*1000);
  end

  % calculate start / end times
  c.tstart = ts_syn(1);
  c.tend = ts_syn(end)+(recsize-1)/c.samplerate;
  
  %%% convert datatype before manipulating data
  if ~isempty(a.convertdatafun) && ~strcmp(func2str(a.convertdatafun), class(c.data)),
    disp(['converting to ' func2str(a.convertdatafun) '...']); % usu. single
    c.data = a.convertdatafun(c.data);
  end
        
  %%% invert sign of data if required
  if invert,
    c.data = -c.data;
  end

  %%% trim time more prettily (by sample now, rather than by record)
  if ~isempty(a.timewin)
    tw = a.timewin;
    if tw(1)<c.tstart
      warning('Requested start time (%.3f s) out of range, using data limits (%.3f s)',...
              tw(1),c.tstart);
      tw(1) = c.tstart;
    end

    if tw(2)>c.tend
      warning('Requested end time (%.3f s) out of range, using data limits (%.3f s)',...
              tw(2),c.tend);
      tw(2) = c.tend;
    end

    c = contwin(c, tw, 'samps_bracket');
  end
  
  %%% update data range
  c = contdatarange(c);
  
  %%% data integrity check
  contcheck(c);
