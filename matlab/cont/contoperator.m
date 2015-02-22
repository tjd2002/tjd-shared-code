function c = contoperator(c1, c2, varargin)
% CONTOPERATOR apply an operator to the data in 2 equal-sized cdats
%
%   cout = contoperator(c1, c2, varargin)
%
% Inputs:
%  c1, c2 - cont structs with equal-sized data
%  'op': operator
%      fnhandle or string of any operator that is valid on the data
%  'name': new name for cdat
%  'infix': string to place between channel names in new channel
%
% Outputs:
%  cout - cont struct with result of operation (time/samplerate from c1)

% Tom Davidson <tjd@stanford.edu> 2003-2010
  
  a = struct('op', [],...
             'name', [],...
             'infix', []);
  
  a = parseArgsLite(varargin,a);

  if isempty(a.infix),
    if ischar(a.op)
      a.infix = ['_' a.op '_'];
    else
        % Use function name if a builtin (func2str won't begin with '@')
        funcstr = func2str(a.op);
        if ~isempty(funcstr) && funcstr(1)~='@'
            a.infix = ['_' funcstr '_'];
        else
          a.infix = '_op_';
        end
    end
  end
  
  % test for equal size
  if ~all(size(c1.data) == size(c2.data))
    error('data must be of same size');
  end
  
  if ~strcmp(class(c1.data), class(c2.data))
    warning('two conts have different data types');
  end
  
  if ischar(a.op),
    a.op = eval(['@' a.op]);
  end

  if ~isa(a.op, 'function_handle')
    error(['''op'' must be (or must evaluate to)_a ' ...
           'function_handle']);
  end
  
  c = c1;
  c.data = a.op(c1.data,c2.data);
  
  c = contdatarange(c);
  
  if isempty(a.name)
    c.name = ['(' c1.name ')' a.infix '(' c2.name ')'];
  else
    c.name = a.name;
  end

  if isempty(c1.chanlabels) && isempty(c2.chanlabels)
    c.chanlabels = {};
  else
    % concatenate chanlabels with operator string
    c.chanlabels = strcat('(', c1.chanlabels, ')', a.infix, '(', c2.chanlabels, ')');
  end
  
  % data integrity check
  contcheck(c);