function c = contcombine(c, cadd, varargin)
% CONTCOMBINE - combine several cont structs that overlap in time into one
%
%  cout = contcombine(c, cadd, [name/value pairs])
%
% Combines cont structs that overlap in time, but may have different
% sampling times, into a single multi-channel cont struct. Resamples and
% interpolates signals so that all channels have same sampling rate and
% sample times.
%
% Default behavior is to resample and interpolate the cont structs in 'cadd'
% so that they match 'c', but a new time basis can also be provided for the
% output cont struct using 'timewin', 'samplerate', and 'nsamps' arguments.
%
% Inputs: (* means required, -> indicates default value)
%   * c - a cont struct to which the others are added
%   * cadd - the cont struct or structs (in an array or cell array) to be 
%       combined with c.
%   'name' - new name for cout
%   'timewin' - requested time range for cout. Actual time range will be the
%       overlap of timewin and the times of all provided structs (default
%       is to use largest overlapping time window across structs).
%   'samplerate'/'nsamps' - requested sampling frequency (Hz) or number
%       of samples for cout. Only one can be provided (default is to
%       use samplerate from c).
%   'interp_method' is the interp1 method that is passed to continterp
%       (->'cubic', 'spline', 'linear', 'nearest', etc)
%   'matchfirst' specifies that the sample times in the first cdat 'c' will be
%       used as-is (new timewin can still be specified). Useful to
%       enforce no interpolation on c. (true,->false)
%
%   Infrequently-used options:
%   'resampbeforeinterpk': (see continterp) controls whether to run an
%       antialiasing filter before resampling. Usually this is a *good
%       idea*, so default is true. Can also be an vector of logicals,
%       with 1 value per cont struct.
%   
% Outputs:
%   cout - the resulting cont struct
%
% Example:
%  Combine and resample several LFP channels that have been imported separately
%
%    cdat_lfp = contcombine(cdat_lfp1, {cdat_lfp2 cdat_lfp3}, 'samplerate', 1500)

%  Tom Davidson <tjd@alum.mit.edu> 2003-2010 

  a = struct('timewin', [],...
             'name', [],...
             'match_first', false,...             
             'nsamps', [],...
             'samplerate', [],...
             'resampbeforeinterpk', true,...
             'interp_method', 'cubic');
  
  a = parseArgsLite(varargin,a);
  
  % validate inputs
  
  if isstruct(cadd), 
    % accept array of structs, but need to use cell array internally as
    % some cdats may have different fields.
    cadd = num2cell(cadd);
  end
  
  if ~isempty(cadd) && ~iscell(cadd)
    error('conts to combine must be in a cell array');
  end

  switch numel(a.resampbeforeinterpk)
   case 1
    resampbeforeinterp = repmat(a.resampbeforeinterpk,numel(cadd)+1,1);
   case numel(cadd)+1
    resampbeforeinterp = a.resampbeforeinterpk;
   otherwise
    error(['must provide either single ''resampbeforeinterpk'' value or one ' ...
           'for each cdat']);
  end

  % data integrity check
  contcheck(c);
  for k = 1:numel(cadd)
    contcheck(cadd{k})
  end
  
  tstart = c.tstart;
  tend = c.tend;
  for k = 1:numel(cadd)
    tstart = max([tstart cadd{k}.tstart]);
    tend = min([tend cadd{k}.tend]);
  end
  
  if isempty(a.timewin),
    timewin = [tstart tend];
  else
    timewin = a.timewin;
  end

  timewin(timewin==-Inf) = tstart;
  timewin(timewin==Inf) = tend;

  
  if a.match_first
    if ~isempty(a.samplerate) || ~isempty(a.nsamps)
      error(['Can''t request match_first and samplerate/nsamps or resampbeforeinterpk(1) (leave ' ...
             'empty to use values from first cdat ''c'')']);
    end
    % don't resample or interpolate, just use requested samples 
    disp('Matching first cdat, no interp or resamp');
    resampbeforeinterp(1) = false;
    c = contwin(c,timewin,'samps_within');
    timewin = [c.tstart c.tend];
  
  else
    if isempty(a.nsamps) && isempty(a.samplerate);
      samplerate = c.samplerate;
    else
      samplerate = a.samplerate;
    end

    if sum([~isempty(samplerate) ~isempty(a.nsamps)]) ~= 1,
      error('exactly one of samplerate/nsamps must be provided');
    end

    % if requested timewin, nsamps or samplerate has changed, interp c
    if all(timewin ~= [c.tstart c.tend]) ||...
          (~isempty(samplerate) && samplerate ~= c.samplerate) ||...
          (~isempty(a.nsamps) && a.nsamps ~= size(c.data,1)),
      c = continterp(c, 'timewin', timewin,...
                     'method', a.interp_method,...
                     'nsamps', a.nsamps,...
                     'resampbeforeinterp', resampbeforeinterp(1),...
                     'samplerate', samplerate);
    end
  end
  
  %%%%%%
  % we now have a reference cdat ('c') with correct timewin, nsamps and samplerate
  %%%%%%
  
  % nsamps in c is what we want to match:
  a.nsamps = size(c.data, 1);
  samplerate = [];
  
  % pre-allocate data array

  disp(sprintf('Pre-allocating memory for combined cdat...', k, length(cadd)));
  dtype = class(c.data);
  newdat = zeros(a.nsamps, numel(cadd)+1, dtype);
  newdat(:,1) = c.data;
  c.data = newdat;
  clear newdat;
  
  for k = 1:length(cadd)

    disp(sprintf('Combining cdat %d of %d...', k+1, numel(cadd)+1));
  
  newname = c.name; 

    if all(timewin ~= [cadd{k}.tstart cadd{k}.tend]) ||...
          (~isempty(samplerate) && samplerate ~= cadd{k}.samplerate) ||...
          (~isempty(a.nsamps) && a.nsamps ~= size(cadd{k}.data,1)),
      cadd{k} = continterp(cadd{k}, 'timewin', timewin,...
                           'nsamps', a.nsamps,...
                           'resampbeforeinterp', resampbeforeinterp(k+1),...
                           'method', a.interp_method);
    end
    
    % concatenate data
    c.data(:,k) = cadd{k}.data;
    
    if ~isempty(c.chanvals) && ~isempty(cadd{k}.chanvals),
      c.chanvals = [c.chanvals cadd{k}.chanvals];
    else
      c.chanvals = [];
    end
    
    if ~isempty(c.chanlabels) && ~isempty(cadd{k}.chanlabels),
      c.chanlabels = [c.chanlabels cadd{k}.chanlabels];
    else
      c.chanlabels = [];
    end
    
    % get extent of data
    c.datarange = vertcat(c.datarange,...
                          cadd{k}.datarange);

    % keep max timestamp error
    c.max_tserr = max(c.max_tserr, cadd{k}.max_tserr);

    % hard to make sense of once we've changed time base for some channels
    c.nbad_start = NaN;
    c.nbad_end = NaN;
    
    % concatenate names
    newname = [newname '&' cadd{k}.name];
  
    % keep nlx_info around if it's present
    try
        c.nlx_info(k+1) = cadd{k}.nlx_info;
    end % ignore caught errors
    
  end
  
  if ~isempty(a.name)
    c.name = a.name;
  else
    c.name = newname;
  end
  
  % data integrity check
  contcheck(c);
