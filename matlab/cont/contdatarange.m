function c = contdatarange(c)
% CONTDATARANGE update datarange field of a contdata struct
% 
% [m x 2] array of upper and lower data ranges for m channels

  nbad_start = c.nbad_start;
  if isnan(nbad_start), 
    nbad_start = 0;
  end

  nbad_end = c.nbad_end;
  if isnan(nbad_end), 
    nbad_end = 0;
  end
  
  % if we trigger an out of memory condition, try some less
  % memory-intensive approaches.
  try
    c.datarange = [min(c.data(nbad_start+1:end-nbad_end,:)); ...
                   max(c.data(nbad_start+1:end-nbad_end,:))]';
  catch
    try 
      c.datarange = [min(c.data); max(c.data)]';
    catch
      % failover to NaNs
      c.datarange = nan(size(c.data,2),2);
    end
  end
  