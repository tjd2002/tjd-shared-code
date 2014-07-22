function [SysMode, recTank, recBlock, BlockStart_unixtime] = TDT_GetStatus(xDA)
%TDT_GETSTATUS Get status of running OpenWorkBench instance

TDT_Connect(xDA);

SysMode = xDA.GetSysMode;

if nargout >1,
  recTank = xDA.GetTankName;
end

if nargout >2,
  xTT = actxserver('TTank.X');
  xTT.ConnectServer('Local','Me');
  xTT.OpenTank(recTank,'R');
  recBlock = xTT.GetHotBlock; % Get block being recorded into
  xTT.SelectBlock(recBlock);
  BlockStart_unixtime = xTT.CurBlockStartTime;
  xTT.CloseTank
  xTT.ReleaseServer
end