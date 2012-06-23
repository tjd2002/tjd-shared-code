function [pks_t pks_idx] = contpeaks(c, varargin)
% CONTPEAKS - find times/indexes of peaks, valleys, or extrema
%
% [pks_t pks_idx] = contpeaks(c, [name/value pairs]
%
% Args:
%  c - cont structure (or any structure with a 'data' field and optionally a
%       'samplerate' field for args in seconds)
%  'type' - 'peaks', 'valleys', 'extrema' (default: peaks)
%
%  'thresh' - threshold in data units
%  'threshfn' - threshold function (defaults to > for peaks, < for
%  valleys, <> for extrema)
%
%  'segs' - only return peaks in specified time segments
%
% Outputs:
%  pks_t -  cell array of peak times, one cell per channel
%  pks_idx - cell array of indexes of peaks, one cell per channel

% Tom Davidson <tjd@stanford.edu> 2003-2010
  

  %%%% input argument parsing/checking

  % data integrity check
  contcheck(c);

  a = struct (...
      'type', 'peaks',...
      'threshfn', [],...
      'thresh',[],...
      'segs',[]);
      
  a = parseArgsLite(varargin, a);

  nchans = size(c.data,2);
  
  switch a.type,
   case 'peaks',
    pkidx = localmax(c.data);
    threshfn_def = @gt;
   case 'valleys'
    pkidx = localmax(-c.data);
    threshfn_def = @lt;
   case 'extrema'
    pkidx = localmax(c.data) | localmax(-c.data);
    threshfn_def = @(x,th) gt(x,th) | lt(x,-th);
   otherwise
    error('bad value for ''type''');
  end
  
  if isempty(a.threshfn)
    a.threshfn = threshfn_def;
  end

  if ~isempty(a.thresh)
    % delete peaks that don't meet threshold
    % (slow, but correct :) )
    pkidx = pkidx & a.threshfn(c.data, a.thresh);
  end
    
  for k = 1:nchans
    pks_idx{k} = find(pkidx(:,k));

    % convert indexes to times
    pks_t{k} = (pks_idx{k}-1)./c.samplerate + c.tstart;
    
    if ~isempty(a.segs)
      goodi = (inseg(a.segs, pks_t{k}));
      pks_t{k} = pks_t{k}(goodi);
      pks_idx{k} = pks_idx{k}(goodi);
    end 
  
  end
  

