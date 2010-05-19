function [c filt] = contfilt(c,varargin)
% CONTFILT filter a contdata structure
  
  a = struct(...
      'filt',[],...
      'filtopt', [],...
      'newname',[],...
      'newchanlabels',[],...
      'autoresample', [],...
      'cache', []);
  
  a = parseArgsLite(varargin,a);
  
  %%% filter, if requested
  if ~xor(isempty(a.filt), isempty(a.filtopt))
    error('exactly one of ''filt'' and ''filtopt'' must be provided');
  end

  if isempty(a.autoresample),
    a.autoresample = true;
  end
  
  %%%
  % If user provides a filter, assume it matches the sample rate of the
  % input signal: filter first, then autoresample if requested. 
  %
  % If user provides a 'filtopt' struct, downsample first for efficiency,
  % then design and apply the filter.
  
  if ~isempty(a.filt),
    
    filt = a.filt;
    
    c = subf_contfilter(c,filt);
    
    if a.autoresample,
      c = subf_autoresamp(c,filt.filtopt);
    end
    
  else
    
    if a.autoresample,
      c = subf_autoresamp(c,a.filtopt);
    end

    % if we are making a filter, set the samplerate from the data
    a.filtopt.Fs = c.samplerate;
    filt = mkfilt('filtopt',a.filtopt, 'cache', a.cache);
    
    c = subf_contfilter(c,filt);
  
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
  
function c = subf_contfilter(c,filt)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter data
  
% can't test whether 'isa' 'filt' object, so:
  try
    get(filt.dfilt,'FilterStructure');
  catch
    error('''dfilt'' field of ''filt''  must be a dfilt structure');
  end
  
  dataclass = class(c.data);
  
  filtlen = length(filt.dfilt.Numerator);

  
  % if samplerate is within 1%, use the filters (i.e. 200Hz would be 204Hz) 
  if abs(log10(filt.filtopt.Fs/c.samplerate)) > log10(1.02),
    error('samplerate of filter doesn''t match data');
  end
  
  % get rid of previous states, if any
  reset(filt.dfilt);
  
  if ~filt.dfilt.islinphase
    % not linear phase, use filtfilt to correct;
    % consider linear phase filters, 1-way filtering is faster
    warning(['filter has non-linear phase distortion, group delay ' ...
             'correction may be incorrect']);
    
    %    warning('dfiltfilt will give sharper responses');
    %    error('dfiltfilt not yet implemented')
    %    c.data = dfiltfilt(a.dfilt,c.data);
    % delta nbad_start and nbad_end will be 2x?
  end

  %  else
    % linear phase (can do one-way filt, then correct for grpdelay)
    if ~mod(filtlen,2),
      warning(['filter length is even, group delay compensation will be ' ...
               'off by 1/2 sample']);
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
    
      
    % correct for group delay, zero-pad end (see nbad, below)
    Gd = filt.dfilt.order/2;
    for k = 1:nchans
      c.data(1:end-Gd,k) = c.data(Gd+1:end,k);
      c.data(end-Gd+1:end,k) = 0;
    end
    
    % indicate which samples are potentially bad
    c.nbad_start = c.nbad_start + filtlen;
    c.nbad_end = c.nbad_end + filtlen;
    
    % cast it back to single, or whatever (necessary?)
    c.data = cast(c.data, dataclass);
    
    % calculate new data range
    c = contdatarange(c);
    
%  end


function c = subf_autoresamp(c,filtopt)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'autoresamp' - resample according to properties of filter, if requested
  
% we do this using the design paramters of the filter, now that we have
% those in the filt structure.
  
% $$$   autothresh = 0.01; % threshold for high end of bandpass

  oversampf = 10; % new Fs factor above 2 * highest freq allowed through 
  
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
% $$$     res_f = hiF*oversampf/c.samplerate;

  if ~strcmp(filtopt.filttype,{'bandpass', 'lowpass'}),
    % don't want to downsample a highpass filter
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

