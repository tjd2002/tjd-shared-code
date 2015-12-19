function slmousedown(fig, ax)
% SLMOUSEDOWN - prepare for zooming on mouseup
%   Note this function works with slmousemove.m and slmouseup.m to handle
%   mouse movement in spikeslpos.m
% $Id$

% get panel containing last object clicked on
panel = ancestor(gco,'uipanel');
if isempty(panel),
  return
end

% get sldata
sld = getappdata(panel,'sldata');

% set up for slmousemove/slmouseup
sld.mousedown = true;
sld.showdist = false;

% click location in pixels from screen bot/left
sld.mousedownpix = get(0,'pointerlocation');

% get click type
sel = get(gcf,'selectiontype');
% on double-click, redo last operation
if strcmp(sel, 'open');
  sel = sld.lastsel;
end
sld.lastsel = sel;

% destroy old zoom boxes
if isfield(sld,'zbarh') && all(ishandle(sld.zbarh)),
  delete(sld.zbarh);
end;
sld.zbarh = [];

% which subax clicked in?
for subax = sld.subaxes,
  % gco is clicked-on object (may be image, plot point...)
  if ismember(gco, findobj(subax))
    % get mousedown location/selection type
    cp = get(subax,'CurrentPoint');
    sld.mousedownx = cp(1,1);
    sld.mousedowny = cp(1,2);
    sld.lastsubax = subax;
    break; % we found our subax
  end
end

% store appdata
setappdata(panel,'sldata',sld);

