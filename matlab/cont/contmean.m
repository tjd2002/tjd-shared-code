function c = contmean(c, varargin)
% CONTMEAN calculate mean of channels in a cont struct
%  c = contmean(c, [name/value pairs])
% 
% Inputs:
%  c - cont struct
%  'name' - New name for struct
%
% Outputs:
%  cout - cont struct with mean values for each channel

% Tom Davidson <tjd@stanford.edu> 2003-2010
  
  a = struct(...
      'name', []);
  
  a = parseArgsLite(varargin,a);
      
  c.data = mean(c.data,2);
  
  c = contdatarange(c);

  c.chanvals = [];
  
  c.chanlabels(2,:) = {'_'};
  c.chanlabels = {[c.chanlabels{:} '_mean']};
  
  if isempty(a.name)
    c.name = [c.name '_mean'];
  else
    c.name = a.name;
  end