function slmovkeypress(fig)
% SLMOVKEYPRESS - keypress callback function for slp2dmovie 
%   
% $Id$
  
% get last key pressed
cc = get(fig,'CurrentCharacter');

a = getappdata(fig,'slmovargs');
d = getappdata(fig,'slmovdata');

switch lower(cc)
  
 case [] % as when shift key is pressed
  return
  
 case ' ' % just play movie, don't return
  replay = true;
  
% $$$  case {'+' '='},
% $$$   sla.maxpfs = sla.maxpfs + 1;
% $$$   
% $$$  case {'-' '_'},
% $$$   if sla.maxpfs > 0
% $$$     sla.maxpfs = sla.maxpfs -1;
% $$$   end
 
 otherwise
  return
  
end

if replay,
  slp2dmovie('argstruct',a);
end
