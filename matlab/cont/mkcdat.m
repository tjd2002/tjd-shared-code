function c = mkcdat(varargin)
%MKCDAT make a new cont struct (constructor)
%
%  Private helper function called by imcont* functions  
%
%  Note name is cdat due to conflict with obj/cont cont objects. This
%  should get resolved, when functionality of cont objects and cont
%  structs is merged. For now it's just cosmetic. :)

% Tom Davidson <tjd@stanford.edu> 2003-2010
  
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
             'max_tserr',[],...
             'mean_tserr',[]);
  
  c = parseArgsLite(varargin,c);
  
  [nsamps nchans] = deal([]); %#ok

  if ~isempty(c.data)
    if ndims(c.data) > 2
      error('continuous data may only be 1-D per channel');
    end
    [nsamps nchans] = size(c.data);
    if nchans>nsamps,
      warning('More samples than channels, is data transposed?');
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
  