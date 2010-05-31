function [bouts th minevpeak minpeak maxevvalley maxvalley extrema extrema_t] = ...
      contbouts(c, varargin)
% CONTBOUTS - find above-threshold periods of a signal in a cont struct
%
%    [bouts thresh minevpeak minpeak maxevvalley maxvalley peakvals peakts] = ...
%        contbouts(cont, [name/value pairs]
%
% Inputs: (* means required, -> indicates default value)
%  * cont - cont structure
%   'timeargunits' units for window, mindur, minevdur, and outputs (->'seconds', 'samples') 
%   'datargunits' units for thresh, minpeak, minevpeak  (->'data', 'stdevs').
%
%  One of the following threshold args is required, unless c.data is a
%  logical array):
%  *'thresh' - threshold in either data units or standard devs above/below mean
%   'thresh_fn' - defaults to '@ge' (>=). Any fn that takes 2 args (data and
%       threshold). Try: @gt, @le, @lt, @eq (>, <=, <, =, respectively) 
%
%  Conditions can be applied to detected events. All are optional and
%  default to [], for no restriction.
%   'window' - acceptable gap between single events (<=) to be merged into
%              a compound event (default [] = don't merge any events)
%   'minevdur' - minimum (>=) valid length of a single event
%   'mindur' - minimum (>=) valid length of compound event
%   'minevpeak' - minimum (>=) peak value for a single event to be valid
%   'minpeak' - minimum (>=) peak value for a compound event to be valid
%   'maxevvalley' - maximum (<=) valley value for a single event to be valid
%   'maxvalley' - maximum (<=) valley value for a compound event to be valid
%
%  Options for how bouts are calculated. Optional.
%   'interp' - When true, estimate threshold crossing times by linear
%       interpolation between sample times (Only valid for 'thresh_fn' of
%       >,>=,<, or <=). When false, return the first and last samples for
%       which thresh_fn returns true (default true).
%   'includeedges' - When true, then bouts start or end at the edge
%       of the data if thresh_fn is true for the first/last sample
%       (default false).
%
% Outputs:
%   bouts - m x 2 array of start & end times 
%   thresh, minevpeak, minpeak, maxevvalley, maxvalley - the values used for
%       these inputs during calculation, in data units (useful if theses were
%       specified in stdevs)
%   extrema, extrema_t - the values and times of the max or min values
%       (depending on thresh_fn) within each bout
%
%
% Example: To find peaks of >= 200, plus the surrounding period when the
% signal is still above 50, and when the total duration of the event is
% at least 0.2 seconds:
%
%  [bursts] = contbouts(cdat, 'thresh', 50, 'minpeak', 200, 'mindur', 0.2)

% Tom Davidson <tjd@alum.mit.edu> 2003-2010 

% TODO
%  - multiple channels, return cell array of bouts?

  %%%% input argument parsing/checking

  a = struct (...
      'timeargunits', 'seconds',...
      'datargunits', 'data', ...
      'thresh_fn', @ge,...
      'thresh',[],...
      'interp',true,...
      'includeedges',false,...
      'minevdur', [],...
      'mindur',[],...
      'window',[],...
      'minevpeak',[],...
      'minpeak',[],...
      'maxevvalley',[],...
      'maxvalley',[]);
      
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
  
  if  ~any(strcmp(func2str(a.thresh_fn), {'ge' 'gt' 'le' 'lt'})),
    if a.interp
      error('interp only supported for >, >=, <, <= (@gt, @ge, @lt, @le)');    
    end
    if nargin > 7,
      % we don't know whether to return peaks or valleys
      error('extrema only supported for >, >=, <, <= (@gt, @ge, @lt, @le)');
    end
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
      maxvalley = datmean - (a.maxvalley * datstd);
      maxevvalley = datmean - (a.maxevvalley * datstd);
     case 'data'
      th = a.thresh;
      minpeak = a.minpeak;
      minevpeak = a.minevpeak;
      maxvalley = a.maxvalley;
      maxevvalley = a.maxevvalley;
    end
    
    % threshold the data: a logical array of where test is met
    th_data = a.thresh_fn(c.data, th);
    
  else
    % data is already a logical
    th = 0.5;
    th_data = c.data;
  end
  
  if ~islogical(th_data)
    warning([mfilename ':ThreshFnError'],...
            ['user-supplied thresh_fn does not return logical array, ' ...
             'converting using ''logical''']); %#ok
    th_data = logical(th_data);
  end,
  
  % crossings from true to false or vice versa
  datxi = find(diff(th_data));

  % calculate interpolation 'factor' (fraction of a sample) to
  % add/subtract from each crossing:
  if a.interp,
    if ~islogical(c.data),
      % have to cast to double to avoid losing precision in time due to
      % imprecise c.data datatype (e.g. int types)
      datxi = datxi + double((c.data(datxi)-th)./(c.data(datxi) - c.data(datxi+1)));
    else
      % for logical inputs, report interpolated crossing time as 1/2-way between
      % samples
      datxi = datxi + 0.5;
    end
  end

  % deal with high start/end of data
  startoffset = 0;
  endoffset = 0;
  
  highstart = false;
  if th_data(1)
    % first xing downward, signal started high
    if ~a.includeedges,
      startoffset = 1; 
      warning([mfilename ':EdgeData'], 'Data begins during bout, ignoring first crossing');
    else
      highstart = true;
      datxi = [1 datxi']';
    end
  end
  
  if th_data(end)
    % last xing upward, signal ends high
    if ~a.includeedges,
      endoffset = 1;
      warning([mfilename ':EdgeData'], 'Data ends during bout, ignoring last crossing');
    else
      datxi = [datxi' numel(th_data)]';
    end
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
    
    % select only bouts with peaks/valleys that meet requirements threshold (use
    % most restrictive of both thresholds, as each event is both an event
    % and a compound event)
    if ~isempty([minpeak minevpeak maxvalley maxevvalley]),
      goodi = subf_pkvaltest(c.data,datxi,...
                             max([minpeak minevpeak]), ...
                             min([maxvalley maxevvalley]));
      datxi = datxi(goodi,:);
    end
    
    % if user just wants raw bouts, we've got them already   
    bouts = datxi;
    
  else
    % we need to calculate durations and gaps
    
    % remove all events less than minevdur/minevpeak before any other processing
    % as if they never happened. poof.
    
    % minevpeak
    if ~isempty([minevpeak maxevvalley]),
      goodevi = subf_pkvaltest(c.data,datxi, minevpeak, maxevvalley);
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
  if ~isempty([minpeak maxvalley]),
    goodevi = subf_pkvaltest(c.data, bouts, minpeak, maxvalley);
    bouts = bouts(goodevi,:);
  end      

  %% get peak/valley times & values during bouts, if requested 
  % (it's easier to do this here at the end than to keep the list of peaks in
  % sync with the changing list of valid bouts)
  if nargout >= 7
    % pre-initialize variables
    nbouts = size(bouts,1);
    extrema = NaN(nbouts,1);
    extrema_t = extrema;

    switch func2str(a.thresh_fn)
     case {'ge' 'gt'}
      ext_fn = @max;
     case {'le' 'lt'}
      ext_fn = @min;
     otherwise
      % should never get here due to arg testing
      error('bad function for extrema testing');
    end
    
    for k = 1:nbouts
      bstarti = floor(bouts(k,1));
      bendi = ceil(bouts(k,2));

      [extrema(k) exti] = ext_fn(c.data(bstarti:bendi));
      extrema_t(k) = bstarti + exti -1;
    end
    
    if strcmp(a.timeargunits, 'seconds');
      extrema_t = c.tstart + ((extrema_t-1) ./ c.samplerate);
    end
  end
  
  % if no interp, convert start times to first sample above th. Do
  % this after finding bouts so that bout durations are correct above
  if ~a.interp,
    bouts(:,1) = bouts(:,1) + 1;

    % except in case of a (valid) high start, which we hack to start at 1
    if bouts(1,1) == 2, bouts(1,1) = 1; end
    
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
  

function goodi = subf_pkvaltest(data,bouts,minpeak,maxvalley)
% find index of bout pairs that have a point above minpeak.
  
  nbouts = size(bouts,1);
  goodi = true(nbouts,1);
  
  for k = 1:nbouts,
    % find(x,1,'first') should be faster than max/min, we only care if any pt is
    % above thresh.
    if ~isempty(minpeak)
      goodi(k) = any(find(data(floor(bouts(k,1)):ceil(bouts(k,2)))>= minpeak,1,'first'));
    end
    if ~isempty(maxvalley)
      goodi(k) = goodi(k) & any(find(data(floor(bouts(k,1)):ceil(bouts(k,2)))<= maxvalley,1,'first'));
    end
  end