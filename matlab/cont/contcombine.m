function c = contcombine(c, cadd, varargin)
% CONTCOMBINE - combine several cont structures, interpolating data
%
% default is to resample and interpolate the cdat structs in 'cadd' so
% that they match 'c', but a new 'timewin'/'nsamps'/'samplerate' can also
% be provided. 
%
% 'interp_method' is passed through to continterp
%
% 'matchfirst' specifies that the samples in the first cdat 'c' will be
% used as-is (new timewin can still be specified) with no interpolation
%
% 'resampbeforeinterpk' controls whether the initial samplerate matching
% interpolation is done before interpolating. When a signal contains
% short stretches of data surrounded by NaNs, this can reduce the edge
% effects caused by interpolation (even with interp methods like
% 'linear'). (Note that for 'nearest' interpolation, resampbeforeinterp
% is always false). Cf. continterp.
  
  a = struct('timewin', [],...
             'nsamps', [],...
             'samplerate', [],...
             'resampbeforeinterpk', true,...
             'match_first', false,...
             'interp_method', 'cubic');
  
  a = parseArgsLite(varargin,a);
  
  %%% interp cdats to be combined so they have the same time basis as 'c'
  if isstruct(cadd), 
    cadd = {cadd};
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

  tstart = c.tstart;
  tend = c.tend;
  for k = 1:length(cadd)
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
  
  for k = 1:length(cadd)
    
    if all(timewin ~= [cadd{k}.tstart cadd{k}.tend]) ||...
          (~isempty(samplerate) && samplerate ~= cadd{k}.samplerate) ||...
          (~isempty(a.nsamps) && a.nsamps ~= size(cadd{k}.data,1)),
      cadd{k} = continterp(cadd{k}, 'timewin', timewin,...
                           'nsamps', a.nsamps,...
                           'resampbeforeinterp', resampbeforeinterp(k+1),...
                           'method', a.interp_method);
    end
    
    % concatenate data
    c.data = [c.data cadd{k}.data];
    
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
    c.name = [c.name '&' cadd{k}.name];
  
    % keep nlx_info around if it's present
    try
        c.nlx_info(k+1) = cadd{k}.nlx_info;
    end % ignore caught errors
    
  end
  
