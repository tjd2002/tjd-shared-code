function c = contmean(c, varargin)
% CONTMEAN return mean of channels in a contstruct
  
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