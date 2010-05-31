function c = contfn(c, varargin)
% CONTFN apply a function to the data in a contdata struct
%
%  cout = contfn(c, [name/value pairs])
%
% Inputs: * = required
%  *c - cont struct
%  *'fn' - function to apply to the data (function name or
%      function_handle)
%  'name' - new name for cout
%  'suffix' - string to append to c.chanlabels (default '_fn');
%
% Ouput:
%  cout - cont struct with processed data
  
% Tom Davidson <tjd@stanford.edu> 2003-2010
  
  a = struct(...
      'fn', [],...
      'name', [],...
      'suffix', '');
  
  a = parseArgsLite(varargin,a);

  [nsamps nchans] = size(c.data);  %#ok
  
  if isempty(a.suffix),
    if ischar(a.fn)
      a.suffix = ['_' a.fn];
    else
      a.suffix = '_fn';
    end
  end
  
  if isempty(a.fn)
    error('Must provide a ''fn'' argument');
  end
  
  if ischar(a.fn),
    a.fn = eval(['@' a.fn]);
  end

  if ~isa(a.fn, 'function_handle')
    error(['''fn'' must be (or must evaluate to)_a ' ...
           'function_handle']);
  end
  
  c.data = a.fn(c.data);
  
  c = contdatarange(c);
  
  if isempty(a.name)
    c.name = [c.name a.suffix];
  else
    c.name = a.name;
  end
  
  for k = 1:nchans
    if ~isempty(c.chanlabels) && ~isempty(c.chanlabels{k})
      c.chanlabels{k} = [c.chanlabels{k} a.suffix];
    end
  end