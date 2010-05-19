function c = mkcdat(varargin)
%MKCDAT make a new cdat struct (constructor)
  
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
  
  c = parseArgsLite(varargin,c);
  
  [nsamps nchans] = deal([]); %#ok

  if ~isempty(c.data)
    if ndims(c.data) > 2
      error('continuous data may only be 1-D per channel');
    end
    [nsamps nchans] = size(c.data);
    if nchans>nsamps,
      warning('more samps than channels, is data transposed?');
    end
  end
  

  if ~isempty(nchans)
    
    % pass through chanvals, if given
    if ~isempty(c.chanvals) && length(c.chanvals) ~= nchans,
      error('''chanvals'' must have as many entries as data columns');
    end
    
    % pass through chanlabels, if given
    if ~isempty(c.chanlabels) && length(c.chanlabels) ~= nchans,
      error('''chanlabels'' must have as many entries as data columns');
    end
  
  end
  