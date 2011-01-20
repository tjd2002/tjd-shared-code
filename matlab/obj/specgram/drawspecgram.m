function hs = drawspecgram(varargin)
% DRAWSPECGRAM plot a specgram object

% todo: 
%   -abs/log(data), 1/f, etc
%   -call draw2d with options
  
  hs = [];
  
  a = struct('sg', [],...
             'opt', [],...
             'ax', [],...
             'xlim', [],...
             'whitebg', false);
  
  
  a = parseArgsLite(varargin,a);
  
  if isempty(a.ax),
    a.ax = gca;
  end

  % massage sg data for display: 1/f correction, take log of strength as
  % does matlab's specgram...
  if ~isempty(a.sg.b),
    a.sg.b = a.sg.b(2:end,:);
    a.sg.f = a.sg.f(2:end);
% $$$     for m = 1:size(a.sg.b,1),
% $$$       a.sg.b(m,:) = a.sg.b(m,:)./a.sg.f(m);
% $$$     end
    a.sg.b = 20*log10(abs(a.sg.b)+eps);
  end

  if ~isempty(a.opt.freqrange),
    yl = a.opt.freqrange;
  else
    if ~isempty(a.sg.f)
      yl = a.sg.f([1 end]);
    else
      yl = [0 1];% else we will just draw text anyway
    end
  end
  
  draw2d('ax', a.ax,...
         'cmap', a.opt.cmap,...
         'cmapwin', a.opt.cmapwin,...
         'datalims', a.opt.datalim,...
         ...
         'data', a.sg.b,...
         'xhistctrs', a.sg.t,...
         'yhistctrs', a.sg.f,...
         ...
         'ylog', a.opt.logfreq,...
         'shading', a.opt.shading,...
         'ylim', yl,...
             ...
         'nodatastring', 'no specgram computed, zoom in or change params');
  
  
  if ~isempty(a.sg.b)
    % draw a scale bar for the analysis window width
    h = scalebar('ax', a.ax,...
                 'barlen', a.sg.t_window,...
                 'bartxt', ['window = ' timestringfmt(a.sg.t_window)],...
                 'color', [1 1 1],...
                 'xl', a.xlim);
    
    hs = [hs; h];
    
    % draw the 'ncycleline' corresponding to the frequency with at least n
    % complete cycles in each analysis window. (Lower freqs are less reliable)
    if ~isempty(a.opt.ncycleline),
      hold on;
      ncyclefreq = a.opt.ncycleline/a.sg.t_window;
      h = plot(a.sg.t([1 end]),...
               [ncyclefreq ncyclefreq],...
               '--',...
               'color', [1 1 1]);
      
      hs = [hs; h];
    end

    % draw a legend for the specgram
    h = contdrawlegend('chanlabels', a.sg.label, ...
                       'dispname', a.opt.dispname,...
                       'plottype', 'image',...
                       'cmaptop', [1 1 1],... % white
                       'ax', a.ax);
    hs = [hs; h];
  end
  
  %%% set xlim (usu according to requested timewin_plot)
  if ~isempty(a.xlim),
    xlim(a.ax, a.xlim),
  end
  
  %%% set requested ylims
  if diff(yl) > 0,
    ylim(a.ax, yl);
  else
    warning('requested datalims are non-increasing; ignoring')
  end
