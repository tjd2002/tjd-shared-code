function [P F] = contpsd(c, varargin)
% CONTPSD compute power spectral density estimate across segs of cont data
%
% [P F] = contpsd(c, [name/value pairs]);
% 
%  WARNING: Experimental code, not well-tested. See pwelch for 
% 
% Inputs:
%  c - cont struct
%  'segs' - restrict analysis to subsets of c.data
%  'method', - method used to calculate psd 
%    -'welch': currently only method available. See help pwelch
%  *'window_t' - time window (in seconds) for calculating PSD
%  'nfft' - number of fft points (see pwelch, etc)
%  'detrend' - whether to detrend the data in each seg (default:
%      'constant'), see help detrend for other options
%
% Outputs:
%  P - power spectral density at frequencies in F
%  F - frequencies at which power is estimated
   
  % data integrity check
  contcheck(c);
  
  a = struct(...
      'segs', [],...
      'method', 'welch',...
      'nfft', [],...
      'window_t', [],...
      'overlap_frac', 0.5,...
      'detrend', 'constant');
  
  a = parseArgsLite(varargin,a);

  switch a.method
   case 'welch',
    %ok
   otherwise
    error('Unrecognized psd method')
  end
  
  if all(size(c.data,2) ~= 1)
    error('Use contchans to select a single channel.');
  end

  % use all data if no segs provided
  if isempty(a.segs),
    a.segs = [c.tstart c.tend];
  end

  % calculate window to use
  if isempty(a.window_t)
    error('Must provide ''window_t''.');
  end
  
  % convert window/noverlap to samples
  window = ceil(a.window_t .* c.samplerate);
  noverlap = ceil(a.window_t.*a.overlap_frac.*c.samplerate);
  
  if noverlap<0 || noverlap>window
    error(['Overlap must be non-negative, and no greater than window ' ...
           'length']);
  end
  
  nsegs = size(a.segs,1);
  
  P = [];
  F = [];
  wtsum = 0;

  for k = 1:nsegs
    c_k = contwin(c,a.segs(k,:));
    dat = c_k.data;

    len = size(dat,1);

    % have to skip if bout shorter than requested window
    if len < window,
      continue
    end

    % detrend if requested
    if ~isempty(a.detrend)
      dat = detrend(dat, a.detrend);
    end

    % compute periodogram
    [P_k F] = pwelch(dat,window,noverlap,a.nfft,c_k.samplerate);

    % weight by # of windows in this seg (we'll re-normalize below)
    wt = floor(len-window./(window-noverlap))+1;
    wtsum = wtsum + wt;
    
    if isempty(P),
      P = P_k .* wt;
    else
      P = P + (P_k .* wt);
    end
  end
  
  % divide by number of windows contributing
  P = P./wtsum;
