function [xc lags_t xc_correct xc_correct_count cdat] = contxcorr(cs, varargin)
% CONTXCORR *EXPERIMENTAL* (see mkcorr!) do an xcorr between 2 cdat chans, 2 cdat structs,
% or acorr on 1 chan/struct. Returns array for several bouts. NaNs treated
% as missing values for 'none_correct' method.
%
% -bouts can be overlapping/unordered.
%
% TODO:
% -deal with short windows? 
% -maxlags in t?

% Tom Davidson <tjd@stanford.edu> 2003-2010

  % data integrity check
  contcheck(cs);

  a = struct(...
      'chans', [],...
      'chanlabels', [],...
      'scaleopt', 'unbiased_coeff',...
      'maxlag_t', [],...
      'bouts',[],...
      'detrend', 'constant',... % default to cross-covariance
      'autoresample', true);
  
  a = parseArgsLite(varargin,a);

  if ~iscell(a.detrend) || length(a.detrend) == 1,
    a.detrend = {a.detrend a.detrend};
  end
  
  % make sure we only have 1 or 2 channels of data, combine them into a
  % single, 1 or 2-channel cdat struct.

  % figure out the largest timewin we can use to calculate xcorrs
  timewin = [max([cs.tstart]) min([cs.tend])];

  if isempty(a.bouts),
    a.bouts = timewin;
  else
   % correct for possible clipping during resampling, combining
   maxsample = 1./min([cs.samplerate]);
   maxbout = max(a.bouts(:))+2*maxsample;
   minbout = min(a.bouts(:))-2*maxsample;
   timewin = [max([cs.tstart, minbout]) ....
              min([cs.tend, maxbout])];
  end
  
  if diff(timewin)<=0, 
    error(['no overlap between cdat inputs (or zero diff between first/last ' ...
           'bout edge)']);
  end

  % hack in 'unbiased_coeff' xcorr method
  if strcmp(a.scaleopt, 'unbiased_coeff');
    unbiased_coeff = true;
    a.scaleopt = 'unbiased';
  else
    unbiased_coeff = false;
  end
  
  if strcmp(a.scaleopt, 'none_correct');
    a.scaleopt = 'none';
    biased_correct = true;
  else
    biased_correct = false;
  end
  
  switch length(cs)
   case 1,
    cdat = contwin(contchans(cs, 'chans', a.chans, 'chanlabels', a.chanlabels),...
                   timewin);
    
    if all(size(cdat.data,2) ~= [1 2])
      error('Exactly 1 or 2 channels must be specified')
    end
   case 2,
    if ~isempty(a.chans) || ~isempty(a.chanlabels),
      error(['If 2 cdats are provided, ''chans'' and ''chanlabels'' may ' ...
             'not  be provided (use contchans on inputs)']);
    end
    
    if size(cs(1).data,2) ~= 1 || size(cs(2).data,2) ~= 1,
      error(['If 2 cdats are provided, they must have single data channels ' ...
             '(use contchans on inputs)']);
    end

    % -align samples by interpolation so that xcorr makes sense
    % -use highest sampling rate for shared sampling rate 
    cdat = contcombine(cs(1), cs(2),...
                       'samplerate', max([cs.samplerate]),...
                       'timewin', timewin);
   otherwise
    error('Only 1 or 2 cdat structs can be cross-correlated');
  end
  
  % convert lags to samples
  maxlags = ceil(a.maxlag_t .* cdat.samplerate);
    
  if any(a.bouts(:)>cdat.tend) || any(a.bouts(:)<cdat.tstart),
    warning('some bouts outside of cdat range, ignoring');
    a.bouts = a.bouts(inseg([cdat.tstart cdat.tend], a.bouts),:);
  end
  
  nbouts = size(a.bouts,1);
  
  % preallocate
  xc = zeros(2*maxlags+1,nbouts);
  if biased_correct,
    xc_correct = xc;
    xc_correct_count = xc;
  else
    xc_correct = [];
    xc_correct_count = [];
  end

  for k = 1:nbouts
    cdat_win = contwin(cdat,a.bouts(k,:));
    A = cdat_win.data(:,1);

    if size(cdat_win.data,2) == 2, % xcorr
      B = cdat_win.data(:,2);

      if ~isempty(a.detrend{1})
        A = detrend(A,a.detrend{1});
      end
      if ~isempty(a.detrend{2})
        B = detrend(B,a.detrend{2});
      end

      if biased_correct
        xc_correct_this = xcorr(~isnan(A),~isnan(B), maxlags, a.scaleopt);
        xc_correct(:,k) = xc_correct(:,k) + xc_correct_this;
        xc_correct_count(:,k) = xc_correct_count(:,k) + xc_correct_this~=0;
        
        
        % NaN treated as missing values: don't increment correction
        % counter, no correlation since value = 0;
        A(isnan(A)) = 0;
        B(isnan(B)) = 0; 
      end

      [xc(:,k) lags] = xcorr(A, B, maxlags, a.scaleopt);

      
    else % acorr

      if ~isempty(a.detrend{1}),
        A = detrend(A,a.detrend{1});
      end
      
      [xc(:,k) lags] = xcorr(A, maxlags, a.scaleopt);
      
      if biased_correct
        this_xc_correct = xcorr(ones(size(A)), maxlags, a.scaleopt);
        xc_correct(:,k) = xc_correct(:,k) + this_xc_correct;
        % keep track of how many distinct bouts contributed to each time lag
        xc_correct_count(:,k) = xc_correct_count(:,k) + this_xc_correct>0;
      end
      
    end
    
    if unbiased_coeff && any(xc(:,k)) % avoid div/0
                                      % scale unbiased xcorr so that central peak is exactly 1
      xc(:,k) = xc(:,k) ./ xc(maxlags+1,k);
    end
    
  end

  if unbiased_coeff && exist('cdat_win', 'var') && size(cdat_win.data,2) > 1
    warning(['correlation values corrected incorrectly for xcorr (set ' ...
             'to 1 for lag=0)']);
  end

  
  % calculate lags in seconds
  lags_t = lags ./ cdat.samplerate;
  