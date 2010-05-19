function [c chans] = contchans(c,varargin)
% CONTCHANS extract/remove channels from contdata struct
  
  a = struct(...
      'chans',[],...
      'chanlabels',[],...
      'remove', false);
  
  a = parseArgsLite(varargin,a);
  
  [nsamps nchans] = size(c.data); %#ok
  
  if ~isempty(a.chans) && ~isempty(a.chanlabels),
    error('can only provide one of chans and chanlabels');
  end
  
  % if chans is empty, use all channels (columns)
  if isempty(a.chans),
    chans = 1:nchans;
  else
    chans = a.chans;
  end
  
  % make a 1x1 cell array out of a string
  if ischar(a.chanlabels) && ~isempty(a.chanlabels),
    a.chanlabels = {a.chanlabels};
  end

  % if chanlabels provided, select them,
  if ~isempty(a.chanlabels),
    chans = [];
    if length(unique(c.chanlabels)) ~= length(c.chanlabels),
      error('repeated chanlabels in cont struct');
    end
    for chanstr = a.chanlabels(:)'
      chans = [chans find(strcmp(chanstr,c.chanlabels))];
    end
  end

  if a.remove, % invert list
    chans = setdiff(1:nchans,chans);
  end
  
  if isempty(chans),
    error('no channels selected to keep');
  end

  % select correct chanvals for data
  if ~isempty(c.chanvals),
    c.chanvals = c.chanvals(chans);
  end
  
  % select correct chanlabels for data
  if ~isempty(c.chanlabels),
    c.chanlabels = c.chanlabels(chans);
  end
  
  % select correct datalims for data
  if size(c.datarange,1)>1
    c.datarange = c.datarange(chans,:);
  end

  % can run out of memory here, believe it or not
  if numel(chans)~=size(c.data,2) || any(chans ~= 1:size(c.data,2))
    c.data = c.data(:,chans);
  end

  % don't change datalims, slow, done in resamp/filt after this call
