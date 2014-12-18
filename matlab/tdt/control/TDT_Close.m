function [xDA] = TDT_dev_init
% Close server connections and release memory used by ActiveX objects. 
% Run before clearing variables.

if ~strcmp(class(xDA),'COM.TDevAcc_X')
  error('xDA is not a TDevAcc.X object.');
end

xDA.CloseConnection; 
xDA.delete; % release memory used by COM/ActiveX object

