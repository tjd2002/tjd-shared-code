function [vals] = contlookup(c, varargin)
% CONTLOOKUP Get values of cont function at particular times (uses interp1/1q)
%
% [vals] = contlookup(c, [name/value pair args]
%
% Inputs:
%  c - cont structure
%  'xt' - times to lookup
%  'method' -  interpolation method. args as for interp1, (default is
%      'nearest', i.e. no interpolation), try 'linear'.
%  'extrap' -  default=Nan , or use any scalar value, or 'extrap' to
%      allow for extrapolation.
%
% Outputs:
%  vals - values at requested times
  
% Tom Davidson <tjd@stanford.edu> 2003-2010
  
%%%% input argument parsing/checking

  a = struct (...
      'xt', [],...
      'method', 'nearest',...
      'extrap', NaN);
 
  a = parseArgsLite(varargin, a);

  if ~isvector(a.xt)
    error('times to lookup must be a single column');
  end
  
  x = (1:size(c.data,1))';
  xi = ((a.xt(:)-c.tstart)*c.samplerate)+1;
  
  switch(a.method)
   case 'linear',
    vals = interp1q(x,c.data,xi);
   otherwise
    vals = interp1(x,c.data,xi,a.method,a.extrap);
  end