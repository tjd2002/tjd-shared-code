function slresize(fig)
% SLRESIZE - redraw spikeslpos panels on resize
% $Id$

% get panel containing last object clicked on
panels = findobj(fig,'type','uipanel');

for panel = panels';
  
  % get slargs
  sla = getappdata(panel,'slargs');

  if isempty(sla) || ...
        strcmp(get(panel,'units'),'pixels'); % weird race condition with
                                             % resize function interrupts
    continue
  end
  
  if isstruct(sla),
    spikeslpos('argstruct', sla, 'panel', panel);
  end
end