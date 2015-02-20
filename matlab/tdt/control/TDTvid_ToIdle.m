function TDT_vididle(xDA, vid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ->Preview (no logging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n\n*** GOING TO IDLE MODE ***'));

dev = xDA.GetDeviceName(0);
if isempty(dev),
    error('Can''t get TDT device name; Try putting workbench in ''Idle'' state, wait a few seconds?');
end

%tdt: TURN OFF VIDACQ
xDA.SetTargetVal([dev '.VidAcq'], 0);

% vid: stop preview and record, if running
pause(0.5) % wait for last frame to arrive.
stoppreview(vid);
pause(0.5);
stop(vid);

if ~isempty(vid.DiskLogger)
  TDT_FrameCnt = xDA.GetTargetVal([dev '.FrameCnt']);
  Vid_FrameCnt = get(vid.DiskLogger, 'FrameCount');
  if TDT_FrameCnt == Vid_FrameCnt,
    disp(sprintf('OK! TDT/Video frame counts match: %d\n', TDT_FrameCnt)); 
  else
    error(sprintf('Frame Count mismatch: TDT = %d; Video file = %d\n', TDT_FrameCnt, Vid_FrameCnt)); 
  end
  % delete logger so that we don't accidentally keep recording into it.
  vid.DiskLogger = [];
else
  warning('No disklogger found--video was not recording.');
end

%tdt: SET MODE TO PREVIEW
TDT_SetMode(xDA, 'idle');

disp(sprintf('   --> DONE\n'));
