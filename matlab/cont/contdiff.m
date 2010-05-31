function c = contdiff(c, varargin)
% CONTDIFF take the diff (derivative) of the signal in a cont struct
% 
%  cout = contdiff(cont, [name/value pairs])
%
% Inputs:
%  'N', take the nth-order derivative (default 1)
%  'name', name for the new cont struct (default: [oldname '_dx' N])
%
% Example:
% 
%    cdat_fr = contdiff(cdat);

% Tom Davidson <tjd@alum.mit.edu> 2003-2010
  
  a = struct(...
      'N', 1,...
      'name', []);
  
  a = parseArgsLite(varargin,a);
      
  [c.data] = diff(c.data, a.N,1);
  c.tstart = c.tstart + (0.5 * a.N/c.samplerate);
  c.tend = c.tend - (0.5 * a.N/c.samplerate);
  
  c = contdatarange(c);
  
  if a.N > 1, 
    derstr = ['_dx' num2str(a.N)];
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