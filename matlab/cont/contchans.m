function [c chans] = contchans(c,varargin)
% CONTCHANS - select channels from cont struct
%
%    [cout chans] = contchans(c, [name/value pairs])
%
% Inputs: (* means required, -> indicates default value)
%   * c - a cont struct (required)
%
%  One of the following two args is required:
%   'chans' - index of channels to select (numeric array)
%   'chanlabels' - labels of channels to select (string or cell array of strings)
%
%  Optional:
%   'remove' - invert the selection; select all channels except those
%        indicated (default false)
%
%
% Outputs:
%  cout - a cont struct containing only the requested channels
%  chans - the indexes of the channels that were selected
%
%
% Example: To select just the channels named 'LFP1' and 'LFP2':
% 
%  cdat = contchans(cdat_lfp, 'chanlabels', {'LFP1' 'LFP2'});

% Tom Davidson <tjd@alum.mit.edu> 2003-2010

  a = struct(...
      'chans',[],...
      'chanlabels',[],...
      'remove', false);
  
  a = parseArgsLite(varargin,a);
  
  [nsamps nchans] = size(c.data); %#ok
  
  if ~isempty(a.chans) && ~isempty(a.chanlabels),
    error('can only provide one of chans and chanlabels');
  end
  
%   % if both args are empty, use all channels (columns)
%   if isempty(a.chans) && isempty(a.chanlabels),
%     chans = 1:nchans;

  % if both args are empty return no channels
  if isempty(a.chans) && isempty(a.chanlabels),
      chans = [];
    
  elseif ~isempty(a.chans),
      chans = a.chans;

  else
      chans = chansfromlabels(c.chanlabels, a.chanlabels);
      
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
