function TDT_dev_connect(xDA)
% Connect to running TDT OpenWorkbench (safe to run multiple times)

% This function opens a connection to a running OpenWorkbench only if
% there is not a currently open one.

if ~strcmp(class(xDA),'COM.TDevAcc_X')
  error('xDA is not a TDevAcc.X object, call TDT_Init to create one.');
end

% When SysMode==0 (Idle state), CheckServerConnection returns 0, even when 
% there is a valid connection, so we attempt to temporarily set SysMode to 
% Standby (1) before checking for a connection.
wasIdle = false;
if xDA.GetSysMode == 0;
  wasIdle = true;
  xDA.SetSysMode(1); % fails silently if no connection, OK
end

if ~xDA.CheckServerConnection;
  out = xDA.ConnectServer('Local');
  if ~out;
    xDA.CloseConnection; % just to be sure?
    error('Can''t connect: Make sure OpenWorkbench is running with a project loaded.');
  end
end

% Return OpenWorkbench to Idle state if we changed it
if wasIdle
  xDA.SetSysMode(0);
end
