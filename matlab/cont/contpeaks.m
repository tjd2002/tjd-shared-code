function [pks_t] = contpeaks(c, varargin)
% CONTPEAKS - find times of peaks, valleys, or zero-crossings
%
% peaks_t = contpeaks(c, [name/value pairs]
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
%  peaks_t -  cell array of peak times, one cell per channel

% Tom Davidson <tjd@stanford.edu> 2003-2010
  

  %%%% input argument parsing/checking

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
    % convert indexes to times
    pks_t{k} = (find(pkidx(:,k))-1)./c.samplerate + c.tstart;
    
    if ~isempty(a.segs)
      pks_t{k} = pks_t{k}(inseg(a.segs, pks_t{k}));
    end 
  
  end
  

