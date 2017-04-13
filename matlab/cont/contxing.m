function [pks_t pks_idx] = contxing(c, varargin)
% CONTXING - find times/indexes of threshold crossings ('xings')
%
% [xing_t xing_idx] = contxing(c, [name/value pairs]
%
% Args:
%  c - cont structure (or any structure with a 'data' field and optionally a
%       'samplerate' field for args in seconds)
%  'type' - 'up', 'down': crossing direction. (Default: up, i.e.
%       from false->true)
%
%  'thresh' - threshold in data units, defaults to 0
%  'equal_xing' - when value == threshold, should we report a crossing?
%      default: true
%
% Outputs:
%  xing_t -  cell array of crossing times, one cell per channel
%  xing_idx - cell array of indexes of crossings, one cell per channel
%
%
% Example:
%
% % get times of upward zero-crossings:
%  xing_t = contxing(c)
%
% % get times of downward crossings of -10
%  xing_t = contxing(c, 'thresh', -10, 'type', 'down')
%

% Tom Davidson <tjd@stanford.edu> 2016

% TODO:
% -add type 'updown' to detect both crossings
% -add 'threshfn' for arbitrary functions of data (but can also use contfn
%  to transform before calling contxing
%  'segs' - only return peaks in specified time segments

  %%%% input argument parsing/checking

  % data integrity check
  contcheck(c);

  a = struct (...
      'type', 'up',...
      'thresh',0,...
      'equal_xing',true);
      
  a = parseArgsLite(varargin, a);

  nchans = size(c.data,2);
  
  switch a.type,
    case 'up',
      if a.equal_xing
        threshfn = @ge;
      else
        threshfn = @gt;
      end
      
    case 'down'
      if a.equal_xing
        threshfn = @le;
      else
        threshfn = @lt;
      end
      
    otherwise
      error('bad value for ''type''');
  end
  
  if isempty(a.thresh) || ~isnumeric(a.thresh)
    error('Must provide non-empty, numeric threshold')
  end

  % find threshold crossings
  dat_xing = diff(threshfn(c.dat),a.thresh)>0;
  
  for k = 1:nchans
    pks_idx{k} = find(dat_xing(:,k));

    % convert indexes to times
    pks_t{k} = (pks_idx{k}-1)./c.samplerate + c.tstart;
    
    if ~isempty(a.segs)
      goodi = (inseg(a.segs, pks_t{k}));
      pks_t{k} = pks_t{k}(goodi);
      pks_idx{k} = pks_idx{k}(goodi);
    end 
  
  end
  

