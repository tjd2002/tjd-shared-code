function [bouts, th] = contboutsbool(c, varargin)
% CONTBOUTSBOOL - find periods in a contdata struct when a given fn is true
%
% [bouts, th] = contboutsbool(cont, [name/value pair args]
%
% Args:
%  fn = anonymous function to test each data sample with
%  cont = contdata structure (or any structure with a 'data' field and
%         optionally a 'samplerate' field for args in seconds)
%
%  'interp' - When true, use estimated crossing time (0.5 samples before
%  up crossing). When false, return first sample above, last sample
%  above.
%
%  'invert' - find bouts when signal is false
%
%  * argunits, minevdur, mindur, and window are as for (and passed
%  through unchanged to) contbouts.
%
%  $Id: contboutsbool.m 2213 2009-08-03 19:38:21Z tjd $

% -algorithm, call contbouts with 'interp' regardless of requested 'interp',
% then round starts up, ends down if interp not requested.
  
%  - multiple channels?
  

  %%%% input argument parsing/checking

  a = struct (...
      'fn', @(data) (data ~= false),... % default to logical(data)
      'argunits', 'seconds',...
      'invert', false,...
      'interp',false,...
      'minevdur', 0,...
      'mindur',0,...
      'window',0);
      
  a = parseArgsLite(varargin, a);

  %% validate args
  
  if isempty(a.fn) || ~isa(a.fn, 'function_handle')
    error('a.fn must be a function_handle');
  end
  
  c.data = a.fn(c.data);
  
  if ~islogical(c.data),
    warning(['contboutsbool function does not return logical array, ' ...
             'converting']);
    c.data = logical(c.data);
  end

  % call contbouts with 'interp', then round if interp not requested
  [bouts th] = contbouts(c, ...
                         'argunits', a.argunits,...
                         'invert', a.invert,...
                         'minevdur', a.minevdur,...
                         'mindur', a.mindur,...
                         'window', a.window,...
                         'interp', true);
                    
  if ~a.interp,
    bouts(:,1) = ceil(bouts(:,1));
    bouts(:,2) = floor(bouts(:,2));
  end
  
  
  
