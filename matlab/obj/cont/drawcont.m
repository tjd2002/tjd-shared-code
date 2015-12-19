function hs = drawcont(varargin)
% DRAWCONT plot a cont object, using prefs from contopt, contdrawopt

  hs = [];
  
  a = struct('cont', [],...
             'opt', [],...
             'plottypeno', [],...
             'whitebg', [],...
             'ax', [],...
             'dpi', [],...
             'xlim', []);
  
  
  a = parseArgsLite(varargin,a);
  
  %%% setup

  if isempty(a.opt),
    a.opt = mkcontdrawopt();
  end

  if isa(a.opt.subsfn, 'function_handle')
    subsfn = a.opt.subsfn;
  else
    switch a.opt.subsfn
     case 'none', 
      subsfn = [];
     case {'max', 'min', 'median'}, 
      subsfn = str2func(a.opt.subsfn);
     case 'extrema'
      subsfn = @subf_extrema;
     otherwise
      error('Unrecognized subsfn: %s', a.opt.subsfn); 
    end
  end
  
  if isempty(a.ax),
    a.ax = gca;
  end

  if isempty(a.whitebg),
      a.whitebg = false;
  end
  
  % get axis size
  oldunits = get(a.ax,'units');
  set(a.ax,'units','pixels');
  axpospix = get(a.ax,'position');
  set(a.ax,'units', oldunits);
  
  % guess at what plot to make if none provided
  if isempty(a.opt.plottype)
    if ndims(a.cont.contdata.data) > 2, 
      error('only 1-D and 2-D data supported');
    else
      a.opt.plottype = 'line';
    end
  end
  
  % set 'linear' (non-'log') scaling
  set(a.ax, 'yscale', 'linear');
  
  % choose a default colororder; rotate it depending on plot number
  if isempty(a.opt.color),
    a.opt.color = hsv(8);
    if ~isempty(a.plottypeno),
      a.opt.color= circshift(a.opt.color,-a.plottypeno+1);
    end
    if a.whitebg,
      % Tone down the value for whitebg
      a.opt.color= rgb2hsv(a.opt.color);
      a.opt.color(:,3) = 0.7;
      a.opt.color= hsv2rgb(a.opt.color);
    end
  end
  
  %%% subsample to speed plotting
  
  if ~isempty(a.dpi)
    set(a.ax, 'units', 'inches');
    axpos_in = get(a.ax, 'position');
    axwidpix = axpos_in(3) * a.dpi;
    set(a.ax, 'units', 'pixels'); % reset
  else
    axwidpix = axpospix(3);
  end
  
  switch a.opt.plottype,
   case {'line' 'area'},
    sampsperpix = 5;
   case 'image',
    sampsperpix = 1;
   case 'stairs',
    sampsperpix = 0.5;
   otherwise,
    error('unrecognized ''plottype''');
  end
  
  nsamps = a.cont.timewini(2) - a.cont.timewini(1);
  
  if nsamps < sampsperpix*axwidpix || ~a.opt.subsample
    % plot all, no subsample
    step = 1;
  else
    % step
    step = floor(nsamps/axwidpix/sampsperpix);
  end
  
  
  % get the time window specified by the timewini indexes
  timewin = a.cont.contdata.tstart + (a.cont.timewini/a.cont.contdata.samplerate);

  % get the data limits to use (ylim in line plot, zlim in image plot)
  if ~isempty(a.opt.datalim),
    datalim = a.opt.datalim;
  else
    ranges = a.cont.contdata.datarange;
    datalim = [min(ranges(:,1)) max(ranges(:,2))];
  end
  
  %%% plot the data
  switch a.opt.plottype,
    
   %% line plot
   case {'line' 'area' 'stairs'}
    
    set(a.ax,'ColorOrder', a.opt.color);
    set(a.ax,'LineStyleOrder', a.opt.linestyle);
    set(a.ax,'NextPlot', 'add'); % preserves colororder on 'plot'
    % plot:
    % data timestamps = cont.tstart + ((indexes-1) * 1/samplerate)
    % vs.
    % data = data(indexes)
    dtimes = a.cont.contdata.tstart + ... 
             (((a.cont.timewini(1):step:a.cont.timewini(2))-1)/ ...
              a.cont.contdata.samplerate);
    
    if ~isempty(subsfn) % block processing

      try
        % blkproc will zero-pad the end, try to get an extra block
        data = a.cont.contdata.data(a.cont.timewini(1): ...
                                    a.cont.timewini(2)+step,:);
        dtimes = [dtimes dtimes(end) + step./a.cont.contdata.samplerate];
      catch
        data = a.cont.contdata.data(a.cont.timewini(1): ...
                                    a.cont.timewini(2),:);
      end

      % do block processing (from image processing toolbox)
      data = blkproc(data,[step,1],  subsfn);

      % blkproc time is properly at center of block
      dtimes = ctrs(dtimes);

      % takes care of ignoring last block, if nec
      data = data(1:numel(dtimes),:);
      
    
    else

      % subsfn = 'none', just choose every nth point
      data = a.cont.contdata.data(a.cont.timewini(1):step:a.cont.timewini(2),:);
    end

    switch(a.opt.plottype)
     
     case 'line',
      h = plot(a.ax, dtimes, data);
     
     case 'area',
      if numel(dtimes)>25000 && strcmp(get(figh, 'renderer'), 'painters'),
        warning(['Large # of data points can choke area with ''painters'' ' ...
                 'renderer']);
        pause(2);
      end

      % by default, area draws 'stacked' area graphs, we want multiple
      % area graphs starting at 0
      for k = 1:size(data,2),
        h(k) = area(a.ax, dtimes, data(:,k));
 
        % area doesn't obey the usual plot color ordering
        ncol = size(a.opt.color,1);
        kmod = mod(k,ncol);
        kmod(kmod==0)=ncol;
        set(h(k), 'facecolor', a.opt.color(kmod,:));
        set(h(k), 'edgecolor', a.opt.color(kmod,:));
      end
     
     case 'stairs'
      % make steps centered on given times
      
      % also need to repeat last value to get stairs to plot it
      sstep = diff(dtimes(1:2));
      sdtimes = [dtimes (dtimes(end)+sstep)] - (sstep/2);

      sdata = data([1:end end],:);
      
      h = stairs(a.ax, sdtimes, sdata);
      
    end
    
% $$$     h = plot(a.ax,...
% $$$              a.cont.contdata.tstart + ... 
% $$$              (((a.cont.timewini(1):step:a.cont.timewini(2))-1)/a.cont.contdata.samplerate),...
% $$$              a.cont.contdata.data(a.cont.timewini(1):step:a.cont.timewini(2),:));
    
    hs = [hs; h];
    
    % contlevels - dashed lines at particular y-values
    % (plot even if no data)
    if ~isempty(a.opt.levels),
      h = plot(a.ax,...
               repmat(timewin(:),1,length(a.opt.levels)),...
               repmat(a.opt.levels,2,1), ...
               '--',...
               'color', a.opt.color(1,:) .* a.opt.levelscolf);

      hs = [hs; h];
    end

    % set requested ylims
    if diff(datalim) > 0,
      ylim(a.ax, datalim);
    else
      %warning('requested datalims are non-increasing; ignoring')
    end
    
    % set requested ydir
    if a.opt.inverty,
      set(a.ax, 'ydir', 'reverse');
    end
    
    %% image plot
   case 'image',
    
    % get image data to plot
    imdata = a.cont.contdata.data(a.cont.timewini(1):step:a.cont.timewini(2),:)';
    
    % y-axis labels
    if ~isempty(a.cont.contdata.chanvals),
      chanvals = a.cont.contdata.chanvals;
    else % just enumerate rows
      chanvals = 1:size(imdata,1);
    end
    
    imh = draw2d('data', imdata,...
               'xhistendedges', timewin,...
               'yhistctrs', chanvals,...
               'ax', a.ax,...
               'cmap', a.opt.cmap,...
               'cmapwin', a.opt.cmapwin,...
               'datalims', datalim,...
               'ydir', 'reverse');
    hs = [hs; imh];

    lims = objbounds(imh);
    ylim(a.ax, lims(3:4));
    
  end

  if a.opt.drawlegend,
    h = contdrawlegend('chanlabels', a.cont.contdata.chanlabels, ...
                       'dispname', a.opt.dispname,...
                       'plottype', a.opt.plottype,...
                       'colororder', a.opt.color,...
                       'whitebg', a.whitebg,...
                       'ax', a.ax);
    
    %'cmaptop', [1 1 1],...
    
    hs = [hs(:); h];
  end
  
  %%% draw scale bar
  if a.opt.drawscalebar,
    if ~isempty(a.opt.scalebarcolor),
      sbcol = a.opt.scalebarcolor;
    else
      sbcol = a.opt.color(1,:);
    end
      
    h = yscalebar('ax', a.ax,...
      'barlen', a.opt.scalebarsize,...
      'bartxt', [num2str(a.opt.scalebarsize) ' ' a.cont.contdata.units],...
      'xl', a.xlim,...
      'whitebg', a.whitebg,...
      'color', sbcol);
    hs = [hs; h];
  end

  
  %%% set xlim (usu according to requested timewin_plot)
  if ~isempty(a.xlim),
    xlim(a.ax, a.xlim),
  end
  
  % remove numbers from x/y axes if requested. Leave ticks
  if ~a.opt.drawxaxis,
    set(a.ax,'xticklabel', []);
  end
  
  if ~a.opt.drawyaxis,
    set(a.ax,'yticklabel', []);
  end
  
function ex = subf_extrema(x)
  
  [ex i] = max(abs(x));
  ex = ex * sign(x(i));
  