function [P F] = contpsd(c, varargin)
% CONTPSD compute power spectral density estimate across segs of cont data

  a = struct(...
      'chans', [],...
      'chanlabels', [],...
      'segs', [],...
      'method', 'welch',...
      'nfft', [],...
      'window_t', [],...
      'detrend', 'constant');
  
  a = parseArgsLite(varargin,a);
  
  % figure out the smallest timewin we can use to calculate xcorrs
  timewin = [max([c.tstart]) min([c.tend])];
  if ~isempty(a.segs),
    timewin = [max(c.tstart, min(a.segs(:))) ....
               min(c.tend, max(a.segs(:)))];
  end
  
  if diff(timewin)<=0, 
    error(['no overlap between cdat inputs (or zero diff between first/last ' ...
           'bout edge)']);
  end
  
  cdat = contwin(contchans(c, 'chans', a.chans, 'chanlabels', a.chanlabels),...
                 timewin);
  
  if all(size(cdat.data,2) ~= 1)
    error('Exactly 1 channel must be specified')
  end
  
  % convert lags to samples
  window = ceil(a.window_t .* cdat.samplerate);
  
  nsegs = size(a.segs,1);
  
  P = [];
  F = [];
  wtsum = 0;

  for k = 1:nsegs
    cdat_k = contwin(cdat,a.segs(k,:));
    dat = cdat_k.data;

    if ~isempty(a.detrend)
      dat = detrend(dat, a.detrend);
    end
      
    len = size(dat,1);

    % have to skip if bout shorter than requested window
    if len < window,
      continue
    end

    [P_k F] = pwelch(dat,window,[],a.nfft,cdat_k.samplerate);

    % not sure about this weight--want it to work same as within pwelch
    wt = floor(len./window);
    wtsum = wtsum + wt;
    
    if isempty(P),
      P = P_k .* wt;
    else
      P = P + (P_k .* wt);
    end
  end
  
  P = P./wtsum;
