function c = contwin(c, timewin, method)
% CONTWIN select a time window from a larger cdat struct
  
  if nargin<3,
    method = 'samps_nearest';
  end
  
  if diff(timewin)<=0
    error('timewindow must be of length > 0');
  end
  
  switch method,
   case 'samps_nearest'
    % need to add 1 since cdat is 1-indexed (i.e. time difference of 0
    % samples means start with sample #1)
    win_samp = round((timewin - c.tstart) * c.samplerate)+1;

   case 'samps_within'
    win_samp(1) = ceil ((timewin(1) - c.tstart) * c.samplerate) +1 ;
    win_samp(2) = floor((timewin(2) - c.tstart) * c.samplerate) +1 ;
    
   case 'samps_bracket'
    win_samp(1) = floor((timewin(1) - c.tstart) * c.samplerate) +1 ;
    win_samp(2) = ceil ((timewin(2) - c.tstart) * c.samplerate) +1 ;
    
   otherwise
    error(['invalid win method: use ''samps_nearest'' (default), ''samps_within'', '...
           ' or ''samps_bracket''.']);
  end
  
  c = contwinsamp(c, win_samp);

  % data integrity check done in contwinsamp
  % contcheck(c);
