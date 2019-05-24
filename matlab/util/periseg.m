function p = periseg(varargin)
% PERISEG Get segment- or event-triggered averages
%
% (cleaned-up version of PERIBOUT)
% TODO: error if try to combine plotshuff with other options like
% choosen, detrend, ...?
  
  a = struct(...
      'ax', [],...
      'segs', [],... (e.g. mua cand events)
      'segs_overall', [],... % segs
      'trigtype', 'segstart',... % or 'segend' or 'event'
      'dur_win', [],...
      'nev_win', [],...
      'nth_ev', [],...
      'events_t', [],... % (e.g. ripples)
      'detrend', '',...
      'segdata_only', false,... % only include data within segs
      'choosen', [],... % choose n random events out of possibles
      'cdat', [],...
      't_pre', [],...
      't_post', [],...
      ...
      'plottype', 'trigave',... % 'trigave', 'trials', 'specgram'
      ...
      'plot_events', false,... %% options for trigave
      'plotci', false,...
      'plotstd', false,...
      'plotsem', false,...
      'plotshuff', false,...
      'tzeroline', true,...
      ...
      'sgopt', [],... %% options for specgram
      'sgdopt', [],...
      ...
      'title', true,...
      'datname', [],... % name of channel (clobbers chanlabel)
      'datalim', [],...
      'linewidth', 2,...
      'color', [1 0 0],...
      'color_ci', [],...
      'meanline', false);

  a = parseArgsLite(varargin,a);

  % save args
  p.args = a;
  p.args.cdat = [];

  if ~isempty(a.detrend) && a.meanline
    error('can''t plot mean when signal has been detrended');
  end
  
  if size(a.cdat.data,2)>1
    error('periseg can only handle single-channel cdats');
  end

  if a.segdata_only
    % to restrict analysis only to within segments, we nan out all other
    % data, then we use nanmean/nansem below to calculate triggered averages
    [~, nonsegsi] = contsegdata(a.cdat, seg_not(a.segs, ...
                                                 'Limits', [a.cdat.tstart a.cdat.tend]));
    for k = 1:size(nonsegsi,1)
      a.cdat.data(nonsegsi(k,1):nonsegsi(k,2),:) = NaN;
    end
  end
  
  if isempty(a.color_ci)
    a.color_ci = max(a.color, [0.7 0.7 0.7]);
  end
  
  p.segs = a.segs;  
  
  % get segs of correct length
  if ~isempty(a.dur_win)
    seglens = diff(p.segs,[],2);
    p.segs = p.segs(seglens > a.dur_win(1) & seglens <= a.dur_win(2),:);
  end
  
  % delete segments that are too close to start/end of data
  nsegs = size(p.segs,1);
  if ~isempty(p.segs)
    p.segs(p.segs(:,1)<(a.cdat.tstart+a.t_pre),:) = [];
    if size(p.segs,1)<nsegs
      warning('Deleted segs from start');
      nsegs = size(p.segs,1);
    end
  end
  
  if ~isempty(p.segs)
    p.segs(p.segs(:,2)>(a.cdat.tend-a.t_post),:) = [];
    if size(p.segs,1)<nsegs
      warning('Deleted segs from end');
      nsegs = size(p.segs,1);
    end
  end
  
  if isempty(p.segs)
      warning('No segs provided to periseg');
      
      if isempty(a.ax)
          figure;
          ax = gca;
      else
          ax = a.ax;
      end
      
      % get axis size
      oldunits = get(ax,'units');
      set(ax,'units','pixels');
      axpospix = get(a.ax,'position');
      set(ax,'units', oldunits);
      
      txtcol = [0.7 0.7 0.7];
      text (axpospix(3)/2, axpospix(4)/2, ... % centered
          'No segments provided',...
          'parent', ax,...
          'units','pixels',...
          'verticalalignment','middle',...
          'horizontalalignment','center',...
          'color', txtcol,...
          'edgecolor', txtcol);
      
      
      return
  end
  
  % get length in samples to calculate
  s_pre = round(a.cdat.samplerate * a.t_pre);
  s_post = round(a.cdat.samplerate * a.t_post);
  nsamps = s_pre + s_post +1;
  
  switch(a.trigtype)
      case 'segstart'
          disp('triggering on seg start');
          trig_string = 'seg start';
          trigs_t = p.segs(:,1);
      case 'segend'
          disp('triggering on seg end');
          trig_string = 'seg end';
          trigs_t = p.segs(:,2);
      case 'event'
          % trigger on events in segments
          [~, nev_per_seg, nth_ev_per_seg] = inseg(p.segs,a.events_t);
          
          % exclude segs outside nev_win
          if ~isempty(a.nev_win)
              p.segs = p.segs(nev_per_seg >= a.nev_win(1) & ...
                  nev_per_seg <= a.nev_win(2),:);
              
      % get new nth_ev with missing segs
      [~, nev_per_seg, nth_ev_per_seg] = inseg(p.segs,a.events_t);
    end
    
    % trigger on nth event
    if isempty(a.nth_ev)
      disp('triggering on all events within segs ');
      trig_string = 'all events within segs';
      trigs_t = a.events_t(nth_ev_per_seg > 0);
    
      if a.plotshuff
        trigs_shuff_t = [];
        for k = 1:size(p.segs,1)
          segk = p.segs(k,:);
          segkdur = diff(segk,[],2);
          trigs_shuff_t = ...
              [trigs_shuff_t ...
               (rand(1, nev_per_seg(k)) * segkdur) + segk(1)];
        end
        trigs_shuff_t = sort(trigs_shuff_t)';
      end
    else
      disp(['triggering on nth event in each segment, (# ' num2str(a.nth_ev) ')']);
      trig_string = sprintf('event # %d within segs',a.nth_ev);
      trigs_t = a.events_t(nth_ev_per_seg == a.nth_ev);
    end
   otherwise        
    error('unrecognized trigtype argument: use start/end');
  end

  if ~isempty(a.choosen)
    if a.choosen > size(trigs_t,1)
      warning('choosen larger than # of triggers, using all');
    else
      disp('using choosen random segs');
      rand_i = randperm(size(trigs_t,1));
      trigs_t = trigs_t(rand_i(1:a.choosen),:);
    end
  end
  
  trigs_samp = round((trigs_t - a.cdat.tstart) * a.cdat.samplerate)+1;
  if a.plotshuff
    trigs_shuff_samp = round((trigs_shuff_t - a.cdat.tstart) * a.cdat.samplerate);
  end
% $$$   % hack off cases where s_pre is 
% $$$   if trigs_samp(1) < s_pre,
% $$$     trigs_samp(1) = [];
% $$$     trigs_t(1) = [];
% $$$   end
% $$$ 
% $$$   if size(a.cdat.data,1)-trigs_samp(end) < s_post,
% $$$     trigs_samp(end) = [];
% $$$     trigs_t(end) = [];
% $$$   end
  
  % initialize
  ntrigs = size(trigs_samp,1);
  p.peritrig_data = zeros(ntrigs, nsamps);
  if a.plotshuff
    p.peritrig_shuff_data = zeros(ntrigs, nsamps);
  else
    p.peritrig_shuff_data = [];
  end
  
  for k = 1:ntrigs
    
    trig = trigs_samp(k);
    p.peritrig_data(k,:) = ...
        a.cdat.data(trig-s_pre:trig+s_post,:);

    if a.plotshuff
      trig_shuff = trigs_shuff_samp(k);
      p.peritrig_shuff_data(k,:) = ...
          a.cdat.data(trig_shuff-s_pre:trig_shuff+s_post,:);
    end
    
    % detrend if requested
    if ~isempty(a.detrend)
        if any(isnan(p.peritrig_data(k,:)))
            if strcmpi(a.detrend , 'constant')
                p.peritrig_data(k,:) = ...
                    p.peritrig_data(k,:)-nanmean(p.peritrig_data(k,:));
            else
                error('Cannot use detrend method other than ''constant'' for data containing nans')
            end
        else % non-NaN data
            p.peritrig_data(k,:) = ...
                detrend(p.peritrig_data(k,:), a.detrend);
        end
    end
    
  end

  p.peritrig = nanmean(p.peritrig_data,1);
  p.peritrig_std = nanstd(p.peritrig_data,0,1); % fails with fieldtrip
  p.peritrig_sem = nansem(p.peritrig_data);
  
  if a.plotshuff
    p.peritrig_shuff = nanmean(p.peritrig_shuff_data);
    p.peritrig_shuff_std = nanstd(p.peritrig_shuff_data,0,1);
    p.peritrig_shuff_sem = nansem(p.peritrig_shuff_data);
  end
  
  p.lags = (-s_pre:s_post) ./ a.cdat.samplerate;

  event_localtimes = [];  
  if a.plot_events
    for k = 1:ntrigs
      timewin = trigs_t(k) + [-a.t_pre a.t_post];
      event_localtimes = [event_localtimes ;...
                         a.events_t(inseg(timewin, a.events_t))...
                         - trigs_t(k)];
    end
  end

  % Get name of channel (for titles)
  dat_name = [];
  if ~isempty(a.datname)
    dat_name = a.datname;
  elseif ~isempty(a.cdat.chanlabels)
    dat_name = a.cdat.chanlabels{1};
  end

  
  %% plot

  % line plots of average signal
  if isempty(a.ax)
    figure; 
    ax = gca;
  else
    ax = a.ax;
  end
  
  switch a.plottype
   case 'trigave'
    % get mean value
    if ~isempty(a.segs_overall)
      p.overall_mean = nanmean(contsegdata(a.cdat,a.segs_overall));
    else
      p.overall_mean = [];
    end
    p.segs_mean = nanmean(contsegdata(a.cdat, p.segs));
    
    % plot 95% confidence interval (+/- ~1.96 * std error of the mean)
    % NOT t-corrected, so shoot me
    plusminus = [];
    if a.plotci
      plusminus = 1.96*p.peritrig_sem;
    elseif a.plotsem
      plusminus = p.peritrig_sem;
    elseif a.plotstd
      plusminus = p.peritrig_std;
    end
    
    if ~isempty(plusminus)
      plusminus(isnan(plusminus)) = 0;
      peritrig_tmp = p.peritrig;
      peritrig_tmp(isnan(peritrig_tmp)) = 0;
      
      patch([p.lags fliplr(p.lags)], ...
            [peritrig_tmp + plusminus,...
             fliplr(peritrig_tmp - plusminus)],...
            a.color_ci, ...
            'edgecolor', 'none',...
            'parent',ax )
    end

    hold(ax, 'on');

    % plot peri-event signal average
    line(p.lags, ...
         p.peritrig, ...
         'parent', ax, ...
         'color', a.color, ...
         'linewidth', a.linewidth);

    
    % plot line at mean signal
    if a.meanline
      if ~isempty(a.segs_overall)
        pmean = p.overall_mean;
      else
        pmean = p.segs_mean;
      end
      
      plot(ax, ...
           [-a.t_pre a.t_post], [pmean pmean],...
           'color', a.color,...
           'linestyle', '--', ...
           'linewidth', a.linewidth);
    end
    
    if isempty(a.datalim)
      if ~isempty(a.cdat.datarange)
%         datalim = a.cdat.datarange;
%       else
        datalim = ylim(ax);
      end
    else
      datalim = a.datalim;
    end
    
    % plot mean with triggers shuffled within events
    if a.plotshuff
      line(p.lags, ...
           p.peritrig_shuff, ...
           'parent', ax, ...
           'color', a.color, ...
           'linewidth', a.linewidth,...
           'linestyle', '--');
    end
    
    % plot line at t = 0
    if a.tzeroline
      plot(ax, [0 0], datalim, 'k--', 'linewidth', a.linewidth);
    end
    
    % title
    if a.title
      title(ax, ['Triggered average of ' dat_name '; n = ' ...
                 num2str(size(p.peritrig_data,1)) newline ,...
                 'trig on ' trig_string], 'interpreter', 'none');
    end

    ylim(ax, datalim);
    
    
   case 'specgram'
     tmp_cdat = imcont('timestamp', p.lags,...
       'data', p.peritrig_data');
    
     wbh = waitbar(0, 'Specgram progress');
    for k = 1:ntrigs
      sg_k = mkspecgram('contdata', tmp_cdat, 'chans',k,...
        'specgramopt', a.sgopt);
      
      % collect data returned(to be averaged)
      sgdat_b(:,:,k) = sg_k.b;
      
      waitbar(k/ntrigs, wbh);

    end
    
    close (wbh);
    
    % copy latest sg result to get times/frequencies right
    sg_mean = sg_k;
    sg_mean.b = mean(sgdat_b,3);
    
    % fudge timestamps to be linear so that draw2d doesn't barf
    sg_mean.t = linspace(sg_mean.t(1), sg_mean.t(end), size(sg_mean.t,2));

    drawspecgram('sg', sg_mean, 'opt', a.sgdopt, 'ax', a.ax);
    
    p.sg = sg_k;
        
    if a.title
      title(ax, ['Average specgram of ' dat_name '; n = ' ...
        num2str(size(p.peritrig_data,1)) newline ,...
        'trig on ' trig_string], 'interpreter', 'none');
    end
    
    
   case 'trials'

       draw2d('ax', ax,...
           'datalims', a.datalim,...
           'data', p.peritrig_data,...
           'xhistendctrs', [-a.t_pre a.t_post],...
           'yhistctrs', 1:size(p.peritrig_data,1),...
           'ydir', 'reverse');
       yticks(ax, 1:size(p.peritrig_data,1));
       ylabel('Trial #');

       % title
       if a.title
           title(ax, ['Selected trials from ' dat_name ';' newline 'n = ' ...
               num2str(size(p.peritrig_data,1))  ,...
               '; trig on ' trig_string], 'interpreter', 'none');
       end

    case 'traces'
        timestamps = linspace(-a.t_pre, a.t_post, size(p.peritrig_data,2));
        plot(ax,timestamps,p.peritrig_data);   
       
  end
  
    if isempty(a.datalim)
      if ~isempty(a.cdat.datarange)
%         datalim = a.cdat.datarange;
%       else
        datalim = ylim(ax);
      end
    else
      datalim = a.datalim;
      end
    
    
    ylim(ax, datalim);
    

  
  % plot times of events
  if a.plot_events
    ee = -a.t_pre:0.01:a.t_post;
    hh = histc(event_localtimes, ee);
    plot(ee,hh./max(hh)./20);
    % scatter(ax, event_localtimes, zeros(length(event_localtimes),1),'>');
  end

  xlabel(ax, 'Time since trigger (s)');
  xlim(ax, [-a.t_pre a.t_post]);

  set(ax, 'layer', 'top');
  