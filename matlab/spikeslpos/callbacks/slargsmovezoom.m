function sla = slargsmovezoom(sla,sld,varargin)
%SLARGSMOVEZOOM - generate new slargs for moves/zooms  
  
% $Id$
  
  args = struct(...
      'movef',0,...
      'movepix',0,...
      'zoomf',1);
  
  args = parseArgsLite(varargin,args);
  
  if args.movef && args.movepix,
    error('Only one of movef and movepix can be specified');
  end
  
  xwid = sld.xend - sld.xstart;
  tmid = mean([sld.xend sld.xstart]);
 
  %%% Move (Pan)
  if args.movef ~= 0,
    args.movepix =  args.movef * xwid;
  end
  
  if args.movepix ~= 0,
    sla.timewin = [(sld.xstart + args.movepix) ...
                   (sld.xend + args.movepix)];
  end
  
  %%% Zoom
  if args.zoomf ~= 1,
    sla.timewin = [tmid - 0.5*(xwid/args.zoomf) ...
                   tmid + 0.5*(xwid/args.zoomf)];
  end
  
  