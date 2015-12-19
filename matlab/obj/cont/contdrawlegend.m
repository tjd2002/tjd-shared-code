function h = contdrawlegend(varargin)
% DRAWCONTLEGEND show cont name
  
  h = [];
  
  a = struct('chanlabels', [],...
             'dispname', [],...
             'plottype', [],...
             'colororder', [],...
             'cmaptop', [],...
             'whitebg', [],...
             'ax', []);

  a = parseArgsLite(varargin,a);

  if ~iscell(a.chanlabels),
    a.chanlabels = {a.chanlabels};
  end
  
  if isempty(a.whitebg),
      whitebg = true;
  end
  
  nchans = size(a.chanlabels,2);
  tstring = [];
 
  % get axis size
  oldunits = get(a.ax,'units');
  set(a.ax,'units','pixels');
  axpospix = get(a.ax,'position');
  set(a.ax,'units', oldunits);
  
  if ischar(a.dispname), % use ischar not isempty so that '' gives no label
    tstring = a.dispname;
    
  elseif ~isempty(a.chanlabels)
    % we have chanlabels we can use
    if nchans == 1, % use chanlabel for 1 channel
      tstring = a.chanlabels;
      
    else
      tstring = []; % don't draw custom text box, use 'legend'/ylabel
      
      switch(a.plottype)
       case {'line' 'stairs' 'area'},
        h = legend(a.ax,...
                   a.chanlabels,...
                   'location','NW'); % upper left = northwest = 'NW'
        set(h,'edgecolor',[0.4 0.4 0.4])
        if a.whitebg,
          set(h,'textcolor','k');
          set(h,'color','w');
        else
          set(h,'textcolor','w');
          set(h,'color','k');
        end
        set(h,'interpreter','none');
        
       case 'image',
        set(a.ax,'ytick', 1:nchans);
        set(a.ax,'yticklabel',...
                  a.chanlabels);
       otherwise
        error('unrecognized contdisp.plottype');
      end
    end
  end
  
  if ~isempty(tstring),

    switch(a.plottype)
     case 'image',
      if ~isempty(a.cmaptop),
        % for images, use color of max value in colormap for text
        txtcol = a.cmaptop(1,:);
      else
        txtcol = [0.7 0.7 0.7];
      end
     case {'line' 'area'}
      % for lines, use first color for text
      txtcol = a.colororder(1,:);
     
     otherwise
      error('unrecognized contdisp.plottype');
    end

    if a.whitebg,
      bgcol = 'w';
% $$$       bgcol = rgb2hsv(txtcol);
% $$$       bgcol(2) = bgcol(2) * 0.1;
% $$$       bgcol = hsv2rgb(bgcol);
    else
      bgcol = txtcol .* 0.2;
    end
    
    h = text(10,axpospix(4)-10, ...
             tstring,...
             'parent', a.ax,...
             'units','pixels',...
             'verticalalignment','top',...
             'color', txtcol,...
             'backgroundcolor',bgcol,...
             'edgecolor', txtcol,...
             'interpreter','none');
  end
  