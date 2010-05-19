function g = getHeaderGains(fh)
% GETGAINS helper fn for mwlIO filehandle to get gains from eegfile
  
  if ~isa(fh,'mwleegfile') && ~isa(fh,'mwlwaveformfile');
    error('Can only get gains from a mwlfixedrecordfile (eegfile, spikefile...');
  end
  
  fhdr = get(fh,'header');

  nchans = str2num(fhdr(2).nchannels);
  
  for k = 1:nchans
    g(k) = str2double(...
        fhdr(2).(['channel ' num2str(k-1) ' ampgain']));
  end
  
  