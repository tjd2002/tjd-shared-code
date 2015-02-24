function c = contdatarange(c)
% CONTDATARANGE private function to update datarange field of a cont struct
%  c = contdatarange(c)

% Tom Davidson <tjd@alum.mit.edu> 2003-2010

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
  
  if isempty(c.datarange) % as when there are no good samples
      c.datarange = nan(size(c.data,2),2);
  end