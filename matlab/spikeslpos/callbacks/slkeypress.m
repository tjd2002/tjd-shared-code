function slkeypress(fig)
% SLKEYPRESS - keypress callback function for spikeslpos 
%   
% $Id$
  
% set small factor to move on shift
sf = 0.2;

% default zoom factor
zf = 2;

% extra 'shift' zoom factor
zfmult = 5;

% # of screens to pan
movef = 0;

% zoom factor
zoomf = 1;

% get last key pressed
cc = get(fig,'CurrentCharacter');


% get panel containing last object clicked on
panels = findobj(gcbf,'type','uipanel');
if isempty(panels),
  return
end

for panel = panels(:)';
  sla = getappdata(panel,'slargs');
  sld = getappdata(panel,'sldata');

  switch lower(cc)

   case [] % as when shift key is pressed
    return
    
    
    %%% Pan

   case ']' % one screen right
    movef = +1;

   case '[' % one screen left
    movef = -1;

   case '}' % smallfactor screens right
    movef = sf;
    
   case '{' % smallfactor screens left
    movef = -sf;

    
    %%% Zoom  
   case '.' % zoom in by zoomfactor
    zoomf = zf;

   case '>' % shift-zoom in by big zoomfactor
    zoomf = zfmult*zf;
    
   case ',' % zoom out by zoomfactor
    zoomf = 1/zf;
    
   case '<' % shift-zoom out by big zoomfactor
    zoomf = 1/(zfmult*zf);

    
    %%% Viewlists  
   case 28 % left arrow
    if ~isempty(sla.viewlist) && sla.viewlisti > 1,
      sla.viewlisti = sla.viewlisti-1;
      % we use viewlist only if no timewin is provided
      sla.timewin = [];
    end

   case 29 % right arrow
    if ~isempty(sla.viewlist) && sla.viewlisti < size(sla.viewlist,1);
      sla.viewlisti = sla.viewlisti+1;
      % we use viewlist only if no timewin is provided
      sla.timewin = [];
    end

    % y zoom toggle
% $$$   case 30 % up arrow
% $$$     sla.posrange = 'fwd';
% $$$ 
% $$$   case 31 % down arrow
% $$$     sla.posrange = 'rev';
% $$$ 
% $$$   case ' ' % space bar
% $$$     sla.posrange = 'all';

    
    %%% Raster Stuff
    
   case {'b', 'B'} % toggle raster Background
    for k = 1:length(sla.plots),
      if isfield(sla.plots{k}, 'rasteropt'),
        sla.plots{k}.rasteropt.rasterbg = ~sla.plots{k}.rasteropt.rasterbg;
      end
    end
% $$$     sla.rasterbg = ~sla.rasterbg;
    
    % plot more or fewer rasters per cell
   case {'+' '='},
    sla.maxpks = sla.maxpks + 1;
    
   case {'-' '_'},
    if sla.maxpks > 0
      sla.maxpks = sla.maxpks -1;
    end
    
    
    %%% Other
    
    % toggle showing/hiding all segments
   case 's'
    nsl = numel(sla.seglists);
    if isempty(sla.plotsegs)||any(sla.plotsegs)
      sla.plotsegs = false(nsl,1);
    else
      sla.plotsegs = true(nsl,1);
    end
    
    % toggle plotting the position in all plots
   case 'p'
    for k = 1:length(sla.plots),
      if isfield(sla.plots{k}, 'parmdrawopt') &&...
            ~strcmp(sla.plots{k}.type, 'parm');
        sla.plots{k}.parmdrawopt.draw = ~sla.plots{k}.parmdrawopt.draw;
      end
    end
    
   case 'h' % full crosshair toggle
    if strcmp(get(fig,'pointer'), 'fullcrosshair'),
      set(fig,'pointer','arrow');
    else
      set(fig,'pointer','fullcrosshair');
    end
    return
    
   case 'r' % randomly permute the cells
    sla.randperm = ~sla.randperm;
    
   case 'f'
    for k = 1:length(sla.plots),
      switch lower(sla.plots(k).display)
       case 'posest',
        sla.plots(k).display = 'raster';
       case 'raster'
        sla.plots(k).display = 'posest';
      end
    end
    
    % play 2-d movie of current plot
   case 'm'
    if ~isfield(sld,'movfigh') || isempty(sld.movfigh)
      sld.movfigh = [];
    end
    
    movmap = sld.ppfmap;
    
    % hack, this should be in slp2dmovie, but will do for now
    if ~isempty(sla.posestscale);
      movmap = movmap .* sla.posestscale;
    end
    
    sld.movfigh = slp2dmovie(...
        'e',evalin('base',sla.ename),...
        'movfigh', sld.movfigh,...
        'pfmaps',movmap,...
        'distbinsize',sla.posestdistbin,...
        'timebinsize',sld.timebin,...
        'timewin',[sld.tstart sld.tend]);
    
    %'pfmaps',cat(3,sld.fpfmap, sld.rpfmap, sld.ppfmap),...

    set(0, 'CurrentFigure', fig);
   otherwise
    return
    
  end

  setappdata(panel,'sldata',sld); % set here in case movie is canceled

  sla = slargsmovezoom(sla,sld,'movef',movef,'zoomf',zoomf);
  
  % redraw
  spikeslpos('argstruct', sla);
end
