function opt = mkcontopt(varargin)
% MKCONTOPT options for creating a cont object for display
%
% $$$    -cont
% $$$      [-filter, resample, etc?... so these can be cached]
  
  opt = struct('filtopt', [],...
               'envopt', [],...
               'autoresample', []);
  
  opt = parseArgsLite(varargin,opt);    

