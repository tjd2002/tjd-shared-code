function ArgStruct=parseArgsLite(args,ArgStruct)
% PARSEARGSLITE Helper function for parsing varargin. 
%
%  ArgStruct= parseArgsLite(args, ArgStruct);
%
% This is my stripped-down version of parseArgs, written by Aslak Grinsted
% and available from the Mathworks file exchange. It doesn't allow for flag
% arguments, aliases, or abbreviations, but it is way faster.
%
% Author: Tom Davidson, tjd@mit.edu
%
%
% Example code:
%
% function out = testfn(varargin)
% defaults = struct('type', 'tree',...
%                   'name', 'oak')
% out = parseArgsLite(varargin, defaults);
% return;
%
%
%  >> testfn('name', 'maple')
%
%  ans = 
%      type: 'tree'
%      name: 'maple'
%
%
% $Id: parseArgsLite.m 2213 2009-08-03 19:38:21Z tjd $

% 'struct', if left to its own devices, will create an *array* of structs when
% given any cell array objects as field values. We have to enclose such
% arguments in a 'value' cell array to avoid this behavior.

if numel(ArgStruct) > 1,
  error(['More than one defaults ArgStruct passed in. Did you use a cell array ' ...
         'with ''struct'' when constructing defaults? If so, use extra braces.']);
end

% Removed Feb 2015: not used in my code, causes problems. (Doesn't error if
% you forget the 'param')
%
% % any non-string arguments at the start of the arglist get put into
% % 'a field called NumericArguments'
% NumArgCount=1;
% while (NumArgCount<=size(args,2))&&(~ischar(args{NumArgCount}))
%     NumArgCount=NumArgCount+1;
% end
% 
% if (NumArgCount>1)
%   NumArgCount=NumArgCount-1;
%   ArgStruct.NumericArguments={args{1:NumArgCount}};
%   args(1:NumArgCount) = [];
% end 


for k = 2:2:length(args), % only applies to 'values'
  if iscell(args{k}),
    args(k) = {args(k)};
  end
end

% parse the new inputs
inargs = struct(args{:});

% overwrite defaults
for fn = fieldnames(inargs)',
  fn = fn{:};
  if ~isfield(ArgStruct, fn),
    error(['Unknown named parameter: ' fn])
  end
  ArgStruct.(fn) = inargs.(fn);
end
