function [c filt] = contfilt(c,varargin)
% CONTFILT - filter data in a cont structure
%
%  [c filt] = contfilt(c, [name/value arg pairs])
%
% Filter all channels in a cont struct using provided filter or filter
% design criteria. By default, downsamples highly oversampled signals before
% filtering for computational efficiency. Compensates for group delay.
%
% Inputs: (* means required, -> indicates default value)
%   * c - a cont struct with data to be filtered
%   *'filtopt' - filter design struct as created by mkfiltopt. See
%      example usage below 
%    'filt' - provide a filt struct as created by mkfilt (or as returned
%       from a previous run of contfilt). Saves filter design time, and
%       useful for maintaining consistency in analysis. 
%    'autoresample' - whether to downsample the signal before
%       low-pass or band-pass filtering. Very conservative and
%       safe--downsampling is to 20x the Nyquist of the end of the high
%       stopband (->true, false).
%
%   Infrequently used inputs:
%    'newname'/'newchanlabels' - by default contfilt will append a string
%       ('_F<filtername>') to the cont.name field and to each channel
%       label. Use these inputs to provide your own naming scheme instead.
%    'cache' - an object cache containing previously designed filters can
%       be searched to save design time. (default []).
%
%   Outputs:
%    cout - cont struct containing filtered (possibly resampled) data
%
%   Example: Filter some data in the theta band:
%
%     fopt = mkfiltopt('filttype', 'bandpass', ...
%                      'F', [4 6 10 12], ... % passband is 6-10 Hz
%                      'name', 'theta')
%
%     cdat_theta = contfilt(cdat, 'filtopt', fopt);

% Tom Davidson <tjd@alum.mit.edu> 2003-2010

  % data integrity check
  contcheck(c);

  a = struct(...
      'filt',[],...
      'filtopt', [],...
      'newname',[],...
      'newchanlabels',[],...
      'autoresample', [],...
      'nonlinphaseok', false,...
      'nodelaycorrect', false,...
      'cache', []);
  
  a = parseArgsLite(varargin,a);
  
  %%% handle input args
  if ~xor(isempty(a.filt), isempty(a.filtopt))
    error('exactly one of ''filt'' and ''filtopt'' must be provided');
  end

  if isempty(a.autoresample),
    a.autoresample = true;
  end
  
  %%%
  % If user provides a filter, assume it matches the sample rate of the
  % input cont, not the downsampled version: filter first, then autoresample
  % if requested.
  %
  % If user provides a 'filtopt' struct, downsample first for efficiency,
  % then design and apply the filter.
  
  if ~isempty(a.filt),
    
    filt = a.filt;
    
    c = subf_contfilter(c,filt,a);
    
    if a.autoresample,
      c = subf_autoresamp(c,filt.filtopt);
    end
    
  else
    
    if a.autoresample,
      c = subf_autoresamp(c,a.filtopt);
    end

    % if we are making a filter, set the samplerate from the data
    a.filtopt.Fs = c.samplerate;

    % see if we already have an appropriate filter in the cache
    filt = mkfilt('filtopt',a.filtopt, 'cache', a.cache);
    
    c = subf_contfilter(c,filt,a);
  
  end
                  
  
  % update name/chanlabels to reflect filtering
  if ~isempty(a.newname),
    c.name = a.newname;
  end
  
  if ~isempty(a.newchanlabels);
    c.chanlabels = a.newchanlabels;
  end
  
  if isempty(a.newname) && isempty(a.newchanlabels),
    c.name = [c.name '_F_' filt.filtopt.name];
    for k = 1:length(c.chanlabels),
      c.chanlabels{k} = [c.chanlabels{k} '_F_' filt.filtopt.name];
    end
  end
  
function c = subf_contfilter(c,filt,a)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter data
  
% can't test whether 'isa' 'dfilt' object, so:
  try
    get(filt.dfilt,'FilterStructure');
  catch
    error('''dfilt'' field of ''filt''  must be a dfilt structure');
  end
  
  dataclass = class(c.data);
  
  filtlen = length(filt.dfilt.Numerator);
  
  % if samplerate is within 1%, use the filters (i.e. 200Hz would be 204Hz) 
  if abs(log10(filt.filtopt.Fs/c.samplerate)) > log10(1.01),
    error(['samplerate of filter doesn''t match data, use a filtopt struct ' ...
           'instead']);
  end
  
  % get rid of previous states, if any
  reset(filt.dfilt);
  
  if ~filt.dfilt.islinphase && ~a.nonlinphaseok,
    
    % not linear phase, could use filtfilt to correct instead of erroring...
    warning(['filter has non-linear phase distortion, delay ' ...
           'correction could be wrong.']); %#ok
    
    %    warning('dfiltfilt will give sharper responses');
    %    error('dfiltfilt not yet implemented')
    %    c.data = dfiltfilt(a.dfilt,c.data);
    % delta nbad_start and nbad_end will be 2x?
    
  end
  
  % linear phase (can do one-way filt, then correct for grpdelay)
  if ~mod(filtlen,2) && ~a.nodelaycorrect,
    warning(['filter length is even, group delay compensation will be ' ...
             'off by 1/2 sample']); %#ok
  end
  
  % if data is 'single', use 'single' arithmetic to save memory
  if strcmp(dataclass,'single');
    filt.dfilt.arithmetic = 'single';
  end
  
  disp('filtering...');
  
  nchans = size(c.data,2);
  
  % call the 'filter' method of dfilt object
  try
    c.data = filt.dfilt.filter(c.data);
  catch
    for k = 1:nchans
      c.data(:,k) = filt.dfilt.filter(c.data(:,k));
    end
  end
  
  
  if ~a.nodelaycorrect
    
    % correct for group delay, zero-pad end (see nbad, below)
    Gd = filt.dfilt.order/2;
    for k = 1:nchans
      c.data(1:end-Gd,k) = c.data(Gd+1:end,k);
      c.data(end-Gd+1:end,k) = 0;
    end
    
    % indicate which samples are potentially bad
    c.nbad_end = c.nbad_end + filtlen;
  end
  
  % bad even if no group delay correction
  c.nbad_start = c.nbad_start + filtlen;
    
%   % cast it back to single, or whatever (necessary?)
%   c.data = cast(c.data, dataclass);
  assert(strcmp(class(c.data),dataclass));
  
  % calculate new data range
  c = contdatarange(c);
  
%  end


function c = subf_autoresamp(c,filtopt)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'autoresamp' - resample according to properties of filter, if requested
  
% we do this using the design paramters of the filter, since we have
% those in the filt structure.
  oversampf = 10; % new Fs factor above 2 * top of highest stop band
  
  %% Old way: actually measure freq response of the filter
% $$$   autothresh = 0.01; % threshold for high end of bandpass
% $$$
% $$$   % get freq response of filter
% $$$   [Fresp Ffreq] = freqz(filt.dfilt, [], c.samplerate);
% $$$   
% $$$   if Fresp(end) > autothresh,
% $$$     warning('Does not appear to be a bandpass/lowpass filter, no resampling');
% $$$   else
% $$$     
% $$$     % find low edge of high stop-band (i.e. highest freq allowed through);
% $$$     hiF = Ffreq(find(Fresp>autothresh,1,'last')+1); 
% $$$     res_f = hiF*2*oversampf/c.samplerate;

  if ~strcmp(filtopt.filttype,{'bandpass', 'lowpass'}),
    % don't want to downsample if using a highpass filter
    return
  end
  
  % resampling factor
  res_f = filtopt.F(end) * 2 * oversampf / c.samplerate;
  
  if res_f >= 1;
    disp('no resampling necessary');
  else        
    % resample
    c = contresamp(c,'resample', res_f);
  end

  % data integrity check
  contcheck(c);
