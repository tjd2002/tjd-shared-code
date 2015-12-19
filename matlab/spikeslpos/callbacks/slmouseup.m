function slmouseup(fig)
% SLMOUSEUP - do zoom for mouse drag zoom
%
% $Id$

%%% constants

% zoom factor (in and out)
zf = 2;

% minimum number of pixels mouse had to move between mousedown and
% mouseup to count as a drag
minxmove = 7;

% number of frames to animate for zf zoom
frames = 5;

% default = redraw plot
replot = true;



% get panel containing last object clicked on
panel = ancestor(gco,'uipanel');
if isempty(panel),
  return
end

% get appdata from figure
sld = getappdata(panel,'sldata');
sla = getappdata(panel,'slargs');

twid = sld.xend - sld.xstart;
tmid = mean([sld.xend sld.xstart]);


if ~isempty(sld.lastsubax),
  % get current pointer location in data coords
  cp = get(sld.lastsubax,'CurrentPoint');
  cpx = cp(1,1);
  cpy = cp(1,2);
  
  % ... and in screen pixel coords
  cppix = get(0,'pointerlocation');

  switch sld.lastsel

   case 'normal' % left button       
    
    % on left-drag, zoom into window
    if abs(cppix(1) - sld.mousedownpix) > minxmove,
      sla.timewin(1) = min(sld.mousedownx, cpx);
      sla.timewin(2) = max(sld.mousedownx, cpx);  
% $$$       
% $$$       % zoom ylims as well if we have dragged in y:
% $$$       if sld.showdist,
% $$$         yl(1) = min(sld.mousedowny, cpy);
% $$$         yl(2) = max(sld.mousedowny, cpy);  
% $$$         for k = 1:length(sla.plots),
% $$$           sla.plots{k}.ylim = yl;
% $$$         end
% $$$       end
      
    % on left-click, zoom in by zf
    else
      sla = slargsmovezoom(sla,sld,'zoomf',zf);
    end

    % delete old zoom bar
    delete(sld.zbarh);
    sld.zbarh = [];

   case 'extend' % middle button
                 % on middle-click, recenter on point
    if abs(cppix(1) - sld.mousedownpix(1)) < minxmove
      sla = slargsmovezoom(sla,sld,'movepix',cp(1)-tmid);
    else % on middle-drag, don't do anything
      replot = false;
    end

    
   case 'alt' % right button
              % on right-click, zoom out by zf
    if abs(cppix(1) - sld.mousedownpix(1)) < minxmove
      sla = slargsmovezoom(sla,sld,'zoomf',1/zf);
    else % on right-drag, don't do anything
      replot = false;
    end
% $$$     
% $$$     % zoom y all the way out on any zoom out
% $$$     for k = 1:length(sla.plots),
% $$$       sla.plots(k).ylim = [];
% $$$     end
% $$$     
  end
end % if sld.lastsubax

% reset mousedrag flag
sld.showdist = false;
sld.mousedown = false;

% save appdata
setappdata(panel,'sldata',sld);

if replot,
  
  if ~sld.prefmousezoom,
    % just redraw, no fancy stuff
    spikeslpos('argstruct', sla);

  else
    % mash xlim to do a pretty zoom
    xl = xlim;

    switch sld.lastsel,

      case 'normal'
        % zoom in, redraw after
        xls = [xl(1):(sla.timewin(1)-xl(1))/frames:sla.timewin(1);
               xl(2):(sla.timewin(2)-xl(2))/frames:sla.timewin(2)];
        xzoom(sla, sld, xls);
        spikeslpos('argstruct', sla);

      otherwise,
        % zoom out/scroll, redraw first
        spikeslpos('argstruct', sla);
        sld = getappdata(panel,'sldata');
        xls = [xl(1):(sld.xstart-xl(1))/frames:sld.xstart; ...
               xl(2):(sld.xend-xl(2))/frames:sld.xend];
        xzoom(sla, sld, xls);

    end
  end % if mousezoom
end % if replot

end % function slmouseup


function xzoom(sla,sld, xls)

% zoom using xlim
set(sld.ax,'xtick',[]);
for j = xls,
  xlim(j');
  for subax = [sld.subaxes sld.ax]
    xlim(subax,j');
  end
  drawnow;
end
set(sld.ax,'xtickmode','auto');
end % subfunction xzoom