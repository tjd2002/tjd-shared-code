function [cout pks pks_t vals vals_t] = contphase(c, varargin)
% CONTPHASE - return a cdat with the wrapped phase of the input
%
% methods:
% 'method's:
%  -interp_peak: interpolate phase between peaks
%  -interp_valley: ditto for valleys 
%  -interp_both; interp peak-valley and valley-peak separately (Klausberger's method)
%  -hilbert (imaginary part of h transform is phase)
  
%
% todo:
%  -new name?
%  -multiple channels (use contlocalmax, loops);
%  -deal with spurious peaks/valleys at start/end of signal (use nbad_*
%  fields?)
  
    
  a = struct(...
      'method', 'hilbert');
  
  a = parseArgsLite(varargin, a);
  
  % init vars that may not be set later
  [pks_idx pks pks_t vals_idx vals vals_t] = deal([]);

  % set up output, NaN out for error-checking
  cout = c;
  cout.data = NaN(size(cout.data));
  
  [nsamps nchans] = size(c.data); %#ok
  
  if nchans > 1, 
    error('contphase currently only supports one channel cdat structs');
  end
  
  disp('calculating phase...');

  switch(a.method)
   case 'hilbert',
    suffix = '_phase_hilb';
    % hilbert seems to be buggy with 'single' data
    cout.data = cast(angle(hilbert(double(c.data))), class(c.data));
    sampwin = [cout.tstart cout.tend];d
   
   case 'interp_peak',
    suffix = '_phase_peaks';
    % localmax works across columns
    % find all minima and maxima
    pks_idx = find(localmax(c.data));

    if length(pks_idx) < 2,
      error('no phase available: need at least 2 peaks in signal');
    end
    
    for k = 1:length(pks_idx)-1
      % linear interpolation between peaks
      cout.data(pks_idx(k):pks_idx(k+1)) = linspace(0,2*pi,diff(pks_idx([k k+1]))+1);
    end

    % we can only calculate phase _between_ peaks
    sampwin = pks_idx([1 end]);

   case 'interp_valley',
    suffix = '_phase_valleys';

    vals_idx = find(localmax(-c.data));
    
    if length(vals_idx) < 2,
      error('no phase available: need at least 2 valleys in signal');
    end
    
    for k = 1:length(vals_idx)-1
      % linear interp between pi-3*pi
      cout.data(vals_idx(k):vals_idx(k+1)) = linspace(pi,3*pi,diff(vals_idx([k k+1]))+1);
    end

    % wrap to 0-2*pi
    hival = cout.data>=2*pi;
    cout.data(hival) = cout.data(hival)-(2*pi);
    
    sampwin = vals_idx([1 end]);
    
   case 'interp_both',
    suffix = '_phase_peakvalley';
    % localmax works across columns
    % find all minima and maxima
    pks_idx = find(localmax(c.data));
    vals_idx = find(localmax(-c.data));

    if length(pks_idx) < 1 || length(vals_idx) < 1,
      error('no phase available: need at least one peak/valley in signal');
    end
          
    
    % the loop below assumes we start and end with a peak.
    % handle other possibilities here, if needed

    if vals_idx(1) < pks_idx(1)
      sampwin(1) = vals_idx(1);
      cout.data(vals_idx(1):pks_idx(1)) = ...
          linspace(pi,2*pi,diff([vals_idx(1) pks_idx(1)])+1);
      vals_idx(1) = [];
    else
      sampwin(1) = pks_idx(1);
    end

    %% ok, now we know are starting with a peak
    
    npks = length(pks_idx);
    nvals = length(vals_idx);
    
    if nvals > npks,
      error('assertion failed: too many valleys');
    end
    
    if nvals == npks,
      sampwin(2) = vals_idx(end);
      cout.data(pks_idx(end):vals_idx(end)) = ...
          linspace(pi,2*pi,diff([pks_idx(end) vals_idx(end)])+1);
    else
      sampwin(2) = pks_idx(end);
    end
    
    %% ok, now we know we are ending with a peak

    %% main loop
    npks = length(pks_idx);
    nvals = length(vals_idx);
    
    for k = 1:npks-1,

      pk = pks_idx(k);
      val = vals_idx(k);
      nextpk = pks_idx(k+1);
      
      cout.data(pk:val) = linspace(0,pi,diff([pk val])+1);
      cout.data(val:nextpk) = linspace(pi,2*pi, diff([val nextpk])+1);
      
    end
    
   otherwise
    error('unsupported phase ''method''');
  end

  % save values, times, if requested
  if nargout > 1,
    
    if isempty(pks_idx);
      pks_idx = find(localmax(c.data));
    end
    
    if isempty(vals_idx);
      vals_idx = find(localmax(-c.data));
    end
    
    % select peaks and valleys that will be in clipped cout
    pks_idx = pks_idx(inseg(sampwin(:)',pks_idx));
    pks = c.data(pks_idx);
    pks_t = cout.tstart + (pks_idx-1)/cout.samplerate;
    
    vals_idx = vals_idx(inseg(sampwin(:)',vals_idx));
    vals = c.data(vals_idx);
    vals_t = cout.tstart + (vals_idx-1)/cout.samplerate;
    
  end

  % hold off on clipping until now, so that we can still use indexes
  cout = contwinsamp(cout, sampwin);
  
  % create new chanlabels

  cout.name = [c.name suffix];
  if ~isempty(c.chanlabels),  
    for k = 1:nchans,
      cout.chanlabels{k} = [c.chanlabels{k} suffix];
    end
  end
  
  % update data range
  cout = contdatarange(cout);
  
  % update units
  cout.units = 'radians';