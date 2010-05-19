function [bouts th minevpeak minpeak peakvals peakts] = contbouts(c, varargin)
% CONTBOUTS - find above-threshold periods of a signal in a contstruct
%
% [bouts, th] = contbouts(cont, [name/value pair args]
%
% Args:
%  cont = cont structure (or any structure with a 'data' field and optionally a
%       'samplerate' field for args in seconds)
%  'timeargunits' (->'seconds'<-, 'samples') units for window, mindur,
%       smoothwin, and outputs
%  'datargunits' (->'stdevs'<-, 'data') units for thresh, minpeak, minevpeak.
%
%  One of the following threshold args is required, unless c.data is a
%  logical array)
%
%  'thresh_fn' - defaults to '@ge' (>=). Any fn that takes 2 args (data and
%      threshold).
%  'thresh' - threshold in either data units or standard devs above/below mean
%
%  'minevdur' - min (>=) length of a single event
%  'mindur' - min (>=) length of compound event
%  'window' - acceptable gap between single events (>=) to make a
%      compound event
%
%  'minevpeak' - minimum peak value of a single event
%  'minpeak' - minimum peak value of a compound event
%
%  'interp' - When true, do linear interpolation of crossing times, When
%      false, return sample before crossing--BEWARE group delay of 0.5
%      smaples: correct yourself or use interp when possible! Only valid
%      for thresh_fn of @gt, @ge, @lt, @le (>,>=,<,<=) (default true)
%
%  $Id: contbouts.m 2213 2009-08-03 19:38:21Z tjd $
%
% -Note: with no mindur/minevdur, bouts can be of zero duration

% todo:
% - optionally include bouts starting and ending at start/end of
%  data. For e.g. runspeed filtering, rather then 'events'
% - with no interp, we have a 1/2 sample delay: we return the samples
%  before the crossing, i.e the time from before the sample goes high to
%  the time after it goes low again. Better estimate of bout duration,
%  but not equiv to >= or >, really. But what's the alternative? when a 
%  signal goes high for one sample, is that a bout of duration 1 or 0?
%  Best is to be explicit, and to interp by default, then let the user
%  round or whatever.
%  - consider rounding effects of converting seconds to samples (warn user)
%  - multiple channels, return cell array of bouts?
%  - return time of peak value during bout  

  %%%% input argument parsing/checking

  a = struct (...
      'timeargunits', 'seconds',...
      'datargunits', 'data', ...
      'thresh_fn', @ge,...
      'thresh',[],...
      'interp',true,...
      'minevdur', [],...
      'mindur',[],...
      'window',[],...
      'minevpeak',[],...
      'minpeak',[]);
      
  a = parseArgsLite(varargin, a);

  if size(c.data,2) > 1,
    error('contbouts only supports one channel at a time, use contchans');
  end
  
  %% validate args
  
  if isempty(a.mindur), a.mindur = 0; end
  if isempty(a.minevdur), a.minevdur = 0; end
  if isempty(a.window), a.window = 0; end
  
  if isempty(a.thresh);
    if ~islogical(c.data)
      error('if data is not logical, a threshold must be provided');
    end
  else
    if islogical(c.data),
      error('can not provide a threshold for logical data');
    end
  end
  
  if a.interp && ~any(strcmp(func2str(a.thresh_fn), {'ge' 'gt' 'le' 'lt'})),
    error('interp only supported for >, >=, <, <= (@gt, @ge, @lt, @le)');    
  end
  
  switch a.timeargunits,
   case 'seconds'
    % convert input params from seconds to samples 
    % don't round, we work in fractional samples in the interp case
    windowsamp = a.window*c.samplerate;
    %    windowsamp_minuseps = windowsamp * (1 - eps(c.samplerate));
    windowsamp_pluseps = windowsamp * (1 + eps(c.samplerate));
    mindursamp = a.mindur*c.samplerate;
    mindursamp_minuseps = a.mindur*c.samplerate * (1 - eps(c.samplerate));
    minevdursamp = a.minevdur*c.samplerate ;
    %    minevdursamp_minuseps = a.minevdur*c.samplerate * (1 - eps(c.samplerate));
    
   case 'samples'
    % args in samples already
    windowsamp = a.window;
    mindursamp = a.mindur;
    minevdursamp = a.minevdur;
    
   otherwise
    error('unrecognized ''timeargunits'' arg. Use ''seconds'' or ''samples''');
    
  end
  
  %%%%%%%%%%%%%%
  % THRESHOLD SIGNAL
  
  if ~islogical(c.data),
    
    % get threshold in samples
    switch a.datargunits
     case 'stdevs',
      datmean = mean(c.data);
      datstd = std(c.data);
      th = datmean + (a.thresh* datstd);
      minpeak = datmean + (a.minpeak * datstd);
      minevpeak = datmean + (a.minevpeak * datstd);
     case 'data'
      th = a.thresh;
      minpeak = a.minpeak;
      minevpeak = a.minevpeak;
    end
    
    % threshold the data: a logical array of where test is met
    th_data = a.thresh_fn(c.data, th);
    
  else
    % data is already a logical
    th = 0.5;
    th_data = c.data;
    
  end
  
  if ~islogical(th_data)
    warning(['contbouts function does not return logical array, ' ...
             'converting']);
    th_data = logical(th_data);
  end,
  
  % crossings from true to false or vice versa
  datxi = find(diff(th_data));
  
  % calculate interpolation 'factor' (fraction of a sample) to
  % add/subtract from each crossing, (only valid for 'ge'???)
  if a.interp,
    if ~islogical(c.data),
      % have to cast to double to avoid losing precision in time due to
      % imprecise c.data datatype
      datxi = datxi + double((c.data(datxi)-th)./(c.data(datxi) - c.data(datxi+1)));
    else
      % for logical inputs, report interpolated crossing time as 1/2-way between
      % samples
      datxi = datxi + 0.5;
    end
  end
  
  % make sure first xing upward, last xing downward 
  startoffset = 0;
  endoffset = 0;
  
  if a.thresh_fn(c.data(1),th),
    % first xing downward, ignore
    startoffset = 1; 
    %    warning([mfilename ':EdgeData'], 'Data begins during bout, ignoring first crossing');
  end
  
  if a.thresh_fn(c.data(end),th),
    % last xing upward, ignore
    endoffset = 1;
    %    warning([mfilename ':EdgeData'], 'Data ends during bout, ignoring last crossing');
  end

  if length(datxi) < (2 + startoffset + endoffset),
    bouts = [];
    return
  end

  % Get rid of unusable xings. Make datxi(1) the first up xing, datxi(end) the
  % last down xing.
  datxi = datxi(1+startoffset:end-endoffset);

  % reshape datxi so that:
  % col 1 = indexes of up xings
  % col 2 = indexes of down xings
  datxi = shiftdim(reshape(datxi,2,[]),1);
  
  if minevdursamp == 0 &&...
        mindursamp == 0 &&... 
        windowsamp == 0,
    
    % select only bouts with peak above requested threshold (use most
    % restrictive of both thresholds, as each event is both an event and
    % a compound event)
    if ~isempty(minpeak) || ~isempty(minevpeak),
      minpki = subf_minpks(c.data,datxi,max([minpeak minevpeak]));
      datxi = datxi(minpki,:);
    end
    
    % if user just wants raw bouts, we've got them already   
    bouts = datxi;
    
  else
    % we need to calculate durations and gaps
    
    % remove all events less than minevdur/minevpeak before any other processing
    % it's as if they never happened. poof.
    
    % minevpeak
    if ~isempty(minevpeak),
      goodevi = subf_minpks(c.data,datxi, minevpeak);
      datxi = datxi(goodevi,:);
    end      
    
    % durgap:
    % col 1 = durations
    % col 2 = gaps
    durgap = subf_getdurgap(datxi);
    
    % minevdursamp
    if minevdursamp > 0,
      goodevi = durgap(:,1) >= minevdursamp;
      datxi = datxi(goodevi,:);
      % recalculate durgap to account for new gaps
      durgap = subf_getdurgap(datxi);  
    end
    

    % if no gap spanning 'window' requested, we don't need the big for loop,    
    if windowsamp == 0;

      % if mindur specified, select those events
      if mindursamp > 0
        goodevi = durgap(:,1) >= mindursamp_minuseps;
        datxi = datxi(goodevi,:);
      end

      bouts = datxi;
      
    else
      % OK, we have to deal with the 'window' and its interaction with 'mindur'
      
      % pre-allocate max size bouts array
      % col 1 = dat index of bout start
      % col 2 = dat index of bout end
      bouts = zeros(size(durgap));

      boutsi = 0;
      foundstart = false;
      
      for k = 1:size(durgap,1),
        
        dur = durgap(k,1);
        gap = durgap(k,2);

        if ~foundstart,
          % looking for a start
          if dur >= mindursamp_minuseps || gap <= windowsamp_pluseps
            % found one
            bout(1) = datxi(k,1); % always start on an up xing
            if gap > windowsamp_pluseps
              % found the end, too
              bout(2) = datxi(k,2); 
              boutsi = boutsi + 1;
              bouts(boutsi,:) = bout;
            else
              % end is at least at the end of the next xing
              foundstart = true;
            end
          end
          
        else 
          % foundstart = true, looking for an end
          if gap > windowsamp_pluseps
            % found the end, now test bout length
            bout(2) = datxi(k,2); % always end on a down xing
            foundstart = false;
            if bout(2) - bout(1) >= mindursamp_minuseps,
              boutsi = boutsi + 1;
              bouts(boutsi,:) = bout;
            else% even grouped, bout not long enough, try again
              bout(:) = 0; 
            end
          end
        end
      end
      
      if foundstart
        % didn't find an end yet, just use last down xing
        bout(2) = datxi(k,2);
        % is that last one any good?
        if bout(2) - bout(1) > mindursamp_pluseps,
          boutsi = boutsi + 1; % keep count
          bouts(boutsi,:) = bout;
        end
      end
      
      % discard empty, preallocated bins
      bouts = bouts(1:boutsi,:);
      
    end

  end   
        
  % if minpeak specified, select only those events
  if ~isempty(minpeak),
    goodevi = subf_minpks(c.data, bouts, minpeak);
    bouts = bouts(goodevi,:);
  end      

  %% get peak times/values during bouts, if requested 
  % (it's easier to do this here at the end than to keep the list of peaks in
  % sync with the changing list of valid bouts)
  peakvals = [];
  peakts = [];
  if nargout >= 5
    nbouts = size(bouts,1);
    peakvals = NaN(nbouts,1);
    peakts = peakvals;

    for k = 1:nbouts
      bstarti = floor(bouts(k,1));
      bendi = ceil(bouts(k,2));

      [peakvals(k) maxi] = max(c.data(bstarti:bendi));
      peakts(k) = bstarti + maxi -1;
    end
    if strcmp(a.timeargunits, 'seconds');
      peakts = c.tstart + (peakts-1) ./ c.samplerate;
    end
  end
  
  % if no interp, convert start crossings to first sample above th. Do
  % this after finding bouts so that bout durations are correct above
  if ~a.interp,
    bouts(:,1) = bouts(:,1) + 1;
  end

  % convert from dat indexes to times if requested
  if strcmp(a.timeargunits, 'seconds');
    % indexes are 1-indexed, tstart is time of 1st sample, not zeroth
    bouts = c.tstart + (bouts-1) ./ c.samplerate;
  end
  
  

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS

function durgap = subf_getdurgap(datxi)
% recompute durgap array
  if isempty(datxi),
    durgap = [];
  else
    % the 'diff' between crossings is alternately a gap or a duration
    % col 1 = samps btw. up xing & dn xing (= durations)
    % col 2 = samps btw. down xing & up xing (= gaps)
    datxi = shiftdim(datxi,1);
    durgap = [diff(datxi(:)); Inf];
    durgap = shiftdim(reshape(durgap,2,[]),1);
  end
  

function minpki = subf_minpks(data,bouts,minpeak)
% find index of bout pairs that have a point above minpeak.
  
  nbouts = size(bouts,1);
  minpki = false(nbouts,1);

  for k = 1:nbouts,
    % find(x,1,'first') should be faster than max, we only care if any pt is
    % above thresh.
    minpki(k) = any(find(data(floor(bouts(k,1)):ceil(bouts(k,2)))> minpeak,1,'first'));
  end