function slmousemove(fig)
% SLMOUSEMOVE - draw zoom bar during mouse movement
%   
% $Id$


% break out of callback if no mousedown (called on mousemove)
% only do anything if in data subax, and on left/middle button drag

% get panel containing last object clicked on
panel = ancestor(gco,'uipanel');
if isempty(panel),
  return
end

sla = getappdata(panel,'slargs');
sld = getappdata(panel,'sldata');

if ~isstruct(sld) || ...
      ~isfield(sld, 'mousedown') || ...
      ~sld.mousedown || ...
      isempty(sld.lastsubax) || ...
      ~any(strcmp(sld.lastsel, {'extend' 'normal' 'alt'})),
  return
end

%% Constants

% threshold for displaying the distance on a drag measure
pixyshowth = 100;

% params for measuring line during middle-drag
draglinecol = [0.25 0.25 0.25];
draglinewid = 6;

%different color highlight/text for zoom/measure
switch sld.lastsel
  case 'normal' % we're zooming
    barcol = [1 1 0]; % yellow
    tlabel = ['\leftrightarrow' sprintf('\n')];
  case {'extend' 'alt'} % we're measuring
    barcol = [0 0 1]; % blue
    tlabel = ['\Delta' sprintf('\n')];
end
barhighlight = (barcol.* 0.05) + get(gcf,'color');


%%% set up vars

% get y limits in axes data units
yl = ylim;

% get current pointer location in subax data coords
cp = get(sld.lastsubax,'CurrentPoint');
cpx = cp(1,1);
cpy = cp(1,2);

% get current pointer location in ax pixels coords
cpf = get(fig,'currentpoint');
axpos = get(sld.ax,'position');
cpa = cpf - axpos(1:2);

% get dist of pointer from mousedown location in screen pixel coords
pixmove = sld.mousedownpix - get(0,'PointerLocation');

% are we doing a y-zoom / y-measure?
% (this is sticky: only unset by slmouseup.)
if abs(pixmove(2)) > pixyshowth;
  sld.showdist = true;
end

% are we dragging L->R or R->L?
if cpx > sld.mousedownx, %drag to left
  rx = sld.mousedownx;
  rw = cpx-sld.mousedownx;
  halign = 'left'; % used for text placement below
elseif cpx < sld.mousedownx %drag to right
  rx = cpx;
  rw = sld.mousedownx-cpx;
  halign = 'right'; % used for text placement below
else % cpx = sld.mousedownx, no rectangle
  return
end

% calc distance dragged in data coords (i.e. cm)
dist = abs(cpy-sld.mousedowny)/100;
vel = dist / rw;


%%% Draw

%draw the zoom bar

dragh = [];

if ~sld.showdist,
  set(fig,'currentaxes',sld.ax);
% $$$ 
% $$$   h = rectangle('position', [rx yl(1) rw yl(2)-yl(1)], ...
% $$$     'facecolor', barhighlight,...
% $$$     'erasemode', 'xor');
% $$$   dragh = [dragh h];

  h = line('xdata', [rx rx], 'ydata', [yl(1) yl(2)], ...
    'color', barcol, ...
    'erasemode', 'xor');
  dragh = [dragh h];

  h = line('xdata', [rx+rw rx+rw], 'ydata', [yl(1) yl(2)], ...
    'color', barcol, ...
    'erasemode', 'xor');
  dragh = [dragh h];

else
  switch sld.lastsel
    case 'normal'
      % we're zooming in x and y, draw a zoombox
      oldgca = get(fig,'currentaxes');
      rectpos = [min([sld.mousedownx cpx]) min([cpy sld.mousedowny]) ...
                 abs(cpx - sld.mousedownx) abs(cpy - sld.mousedowny)];
      for subax = sld.subaxes,
        set(fig,'currentaxes',subax);
        h = rectangle('position', rectpos, ...
          'facecolor', 'k',...
          'edgecolor', barcol, ...
          'erasemode', 'xor');
        dragh = [dragh h];
        set(fig,'currentaxes',oldgca);
      end

    case {'extend' 'alt'}
      % draw a diagonal line in each plot
      oldgca = get(fig,'currentaxes');
      for k = 1:length(sla.plots),
        % test if it's got position as the y-axis
        if any(strcmp(sla.plots{k}.type,{'parmest','parm'})) ||...
              (strcmp(sla.plots{k}.type, 'raster') && ~sla.plots{k}.rasteropt.ranked),
          subax = sld.subaxes(k);
          set(fig,'currentaxes',subax);
          h = line([sld.mousedownx; cpx],[sld.mousedowny; cpy], ...
                   'color', draglinecol,...
                   'erasemode', 'xor',...
                   'linewidth', draglinewid);
          dragh = [dragh h];
        end
      end
      set(fig,'currentaxes',oldgca);
  end  
end

% calculate offset for text box next to pointer
if strcmp(halign, 'right'),
  cpa(1) = cpa(1) - 5;
else
  cpa(1) = cpa(1) + 15;
end
cpa(2) = cpa(2) - 25;

% draw text box next to pointer

[tstring secstring] = timestringfmt(rw);

freqstring = '';
if rw < 0.5
  freqstring = [num2str(1/rw, '%.f') ' Hz'];
end

diststring = {};
velstring = '';

if sld.showdist,
  % add distance info to the textbox
  diststring = {[num2str(dist, '%.2f') ' m'] ...
    ['(' num2str(abs(sld.mousedowny), '%.f') ' - ' num2str(abs(cpy), '%.f') ')']};
  if vel > 2;
    velstring = [num2str(vel, '%.1f') ' m/s'];
  else
    velstring = [num2str(vel, '%.2f') ' m/s'];
  end

end


tstring = {[tlabel tstring secstring] ...
           freqstring ...
           diststring{:} ...
           velstring};

% remove empty lines
tstring = tstring(~strcmp(tstring, ''));

textbgcol = get(gcf,'color');
textcol = 1-textbgcol;

texth = text('string', tstring, ... 
  'fontunits', 'pixels', ...
  'fontsize', 12, ...
  'units', 'pixels', ...
  'position', cpa, ...
  'margin', 4, ...
  'backgroundcolor',textbgcol,...
  'color', textcol, ...
  'erasemode', 'xor');



% consolidate to one handle array
newzbarh = [dragh texth];

% delete the old zoom bar
delete (sld.zbarh);

% store the new zoom bar handles in fig appdata
sld.zbarh = newzbarh;
setappdata(panel,'sldata',sld);

end
