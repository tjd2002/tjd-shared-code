function [SysMode, recTank, recBlock, BlockStart_unixtime] = TDT_SetMode(xDA, modestr)
%TDT_SETMODE Change OpenWorkbench modes

switch lower(modestr)
  case {'r' 'record'};
    mode = 3;
  case {'p' 'preview'};
    mode = 2;
  case {'s' 'standby'};
    mode = 1;
  case {'i' 'idle'};
    mode = 0;
  otherwise
    error('Unrecognized mode: ''%s''. Valid modes are ''record''/''preview''/''standby''/''idle''.',...
      modestr);
end

TDT_Connect(xDA);

out = xDA.SetSysMode(mode);
if ~out
  error('Unable to switch to mode ''%s'' in TDT OpenWorkbench. Check ''Messages'' tab?', modestr);
end

t = tic;
while xDA.GetSysMode ~= mode
  pause(0.1);
  if toc(t) > 5, % wait up to 5 seconds 
    errstr = 'Timeout waiting for mode switch in TDT OpenWorkbench.';
    if mode == 3; % Can't record if blockname prompt is on.
      errstr = [errstr ' (Blockname prompting on?)'];
    end
    error(errstr);
  end
end
%toc(t)
[SysMode, recTank, recBlock, BlockStart_unixtime] = TDT_GetStatus(xDA);

