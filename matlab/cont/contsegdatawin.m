function [data] = contsegdatawin(cont, t, offset)
% CONTSEGDATA get windows of data around time points, return as a matrix
  
  if size(cont.data,2) > 1, 
    error ('only one-channel supported right now, use contchans');
  end
  
  t_samp = round((t - cont.tstart) * cont.samplerate)+1;

  offset_samp = round(offset * cont.samplerate);
  
  % omit time windows that are off edge of available data
  t_samp(t_samp <= -offset_samp(1) | t_samp >= size(cont.data,1)-offset_samp(2)) ...
      = [];
  
  for k = 1:size(t_samp,1),
    
    data(k,:) = cont.data(t_samp(k)+offset_samp(1):t_samp(k)+offset_samp(2));

  end
  
