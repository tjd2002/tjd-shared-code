function newc = contwinsamp(c, win_samp)
% CONTWINSAMP select a window of samples from a larger cdat struct
  
  if diff(win_samp)<0
    error('sample window may not be negative length');
  end
  
  if win_samp(1) < 1,
    error('timewindow start time out of range');
  end
  
  if win_samp(2) > size(c.data,1),
    error('timewindow end time out of range');
  end
  
  newc = c;
  newc.data = c.data(win_samp(1):win_samp(2),:);
  
  % subtract one sample, since we want to add the difference between the
  % previous start (at sample #1) from the current start.
  newc.tstart = c.tstart + ((win_samp(1)-1) / c.samplerate);

  newc.tend = newc.tstart + (diff(win_samp)/ c.samplerate);

  % # of bad samples at edges due to filtering, etc.
  newc.nbad_start = c.nbad_start - win_samp(1);
  if newc.nbad_start <0, 
    newc.nbad_start = 0; 
  end
  
  newc.nbad_end = win_samp(2) - (size(c.data,1)-c.nbad_end) ;
  if newc.nbad_end <0,
    newc.nbad_end = 0;
  end
    
  % calc new data range
  newc = contdatarange(newc);
  
  % max_tserr can't be recalculated, but is an upper bound, so we'll keep
  % it
  
  % data integrity check
  contcheck(c);
