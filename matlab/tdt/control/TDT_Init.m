function [xDA] = TDT_dev_init
% Create a TDevAcc object and connect to running OpenWorkbench

% Use actxserver (COM server object) instead of actxcontrol (activeX
% control, which is a type of COM server) to avoid creating spurious figuer
% windows. For ActiveX OCX files like TDT distributes, Matlab creates 
% in-process servers, so these should be just as fast.
xDA = actxserver('TDevAcc.X');

out = xDA.ConnectServer('Local');
if ~out;
  xDA.CloseConnection; % just to be sure?
  error('Can''t connect: Make sure OpenProject is running with a project loaded.');
end

% report something fun
disp(xDA.GetTankName);



