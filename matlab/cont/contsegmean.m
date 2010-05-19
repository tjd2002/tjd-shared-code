function [stat] = contsegmean(c, varargin)
% CONTSEGMEAN measure mean/std/s.e.m. in windows during segments
%
% Since the neighboring samples in biological time series are often
% correlated (that is, not independent), it is difficult to get an estimate
% of the variance of the mean. 
% 
% We overcome this by only making measurements at intervals that are
% expected to be independent. 
  
  a = struct(...
      'segments', [],...
      'meas_dur', [],...
      'meas_count', [],...
      'meas_interval', []);
  
  a = parseArgsLite(varargin,a);
  
  [nsamps nchans] = size(c.data);  %#ok
  
  switch sum([~isempty(a.meas_count) ~isempty(a.meas_interval)])
   case 0

    % we will use the central 'meas_dur' in each segment
    if isempty(a.segments),
      error('Must provide ''segments'' list');
    end

    if size(a.segments,2) ~= 2,
      error('''segments'' must be an m x 2 array of times');
    end

   case 1,
    % ok!
    
   case 2,
    
    error(['only one of ''meas_count'' and ''meas_interval'' can be ' ...
           'provided']);
  
  end
  
  % select segments of at least meas_dur
  segdurs = diff(a.segments,[],2);
  goodsegi = segdurs >= a.meas_dur;
  if ~all(goodsegi)
% $$$     warning(['Only segments at least as long as ''meas_dur'' will be included ' ...
% $$$              'in calculation of mean']);
    a.segments = a.segments(goodsegi,:);
    segdurs = segdurs(goodsegi);
  end
  
  if ~isempty(a.meas_interval),
    % make the average interval between randomly selected measurements
    a.meas_count = round(sum(segdurs) ./ a.meas_interval);
  end
  
  if ~isempty(a.meas_count),
    % we need to randomly choose windows from the available data. We do this
    % using the inverse CDF method. (To avoid biasing selection towards
    % shorter segments)
    
    CDFi = [0; cumsum(segdurs-a.meas_dur)]; % subtract dur, since we can't start any
                                      % closer than that to the end of each seg.
    CDFmax = CDFi(end);
    randt = rand(a.meas_count,1)*CDFmax;
    
    % map the randt's (which are valid on the segments) back into cont
    % timestamp time
    segstarts = zeros(a.meas_count,1);
    for k = 1:a.meas_count,
      segi = find(randt(k)<CDFi,1,'first')-1;
      segstarts(k) = a.segments(segi,1) + (randt(k)-CDFi(segi));
    end
    
    nsegs = a.meas_count;
  else
    
    % get start times of windows of length dur centered on middle of segments:
    segstarts = mean(a.segments,2) -(a.meas_dur/2);
  
    nsegs = size(segstarts,1);
  end
  
  % convert to samples
  segstarts_samps = round((segstarts - c.tstart) * c.samplerate) + 1;
  dur_samps = round(a.meas_dur * c.samplerate);
  
  segmeans = NaN(nsegs, nchans);

  for k = 1:nsegs;
    segsamps = segstarts_samps(k):(segstarts_samps(k)+dur_samps);
    data = c.data(segsamps,:);
    segmeans(k,:) = mean(c.data(segsamps,:));

    % force mean to be NaN if any data are NaN
    segmeans(k,any(isnan(data))) = NaN;
    
  end
  
  % exclude segments containing NaNs/Infs, with warning.
  if any(isnan(segmeans(:))),
    warning(['segments with any NaN values will not be included in calculation ' ...
             'of mean']);
  end
  
  for j = 1:nchans,
    % don't include NaNs (short segs, segs with NaNs...)
    means = segmeans(~isnan(segmeans(:,j)),j);

    stat.segmeans(:,j) = means;
    stat.mean(j) = mean(means);
    stat.std(j) = std(means);
    stat.se(j) = sem(means);
    stat.n(j) = size(means,1);
  end