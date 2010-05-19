function c = contdiff(c, varargin)
% CONTDIFF take the diff of the signal
  
  a = struct(...
      'N', 1,...
      'name', []);
  
  a = parseArgsLite(varargin,a);
      
  [c.data] = diff(c.data, a.N,1);
  c.tstart = c.tstart + (0.5 * a.N/c.samplerate);
  c.tend = c.tend - (0.5 * a.N/c.samplerate);
  
  c = contdatarange(c);
  
  if a.N > 1, 
    derstr = ['_d' num2str(a.N) 'x'];
  else
    derstr = '_dx';
  end
  
  if isempty(a.name)
    c.name = [c.name derstr];
  else
    c.name = a.name;
  end
  
  for k = 1:length(c.chanlabels),
    c.chanlabels{k} = [c.chanlabels{k} derstr];
  end