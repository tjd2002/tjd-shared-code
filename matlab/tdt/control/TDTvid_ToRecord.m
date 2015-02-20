function TDT_vidrecord(xDA, vid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ->Record into new block, new video file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n\n*** GOING TO RECORD MODE ***'));

src = getselectedsource(vid);
dev = xDA.GetDeviceName(0);
if isempty(dev),
    error('Can''t get TDT device name; Try putting workbench in ''Idle'' state, wait a few seconds?');
end

% %tdt: Save VIDACQ state
% VidOn = xDA.GetTargetVal([dev '.VidAcq']);

% TURN OFF VIDACQ
xDA.SetTargetVal([dev '.VidAcq'], 0);

% vid: stop preview and record, if running
pause(0.5);
stoppreview(vid);
stop(vid);

if ~isempty(vid.DiskLogger)
  TDT_FrameCnt = xDA.GetTargetVal([dev '.FrameCnt']);
  Vid_FrameCnt = get(vid.DiskLogger, 'FrameCount');
  if TDT_FrameCnt == Vid_FrameCnt,
    disp(sprintf('OK! TDT/Video frame counts match: %d\n', TDT_FrameCnt)); 
  else
    warning(sprintf('Frame Count mismatch: TDT = %d; Video file = %d\n', TDT_FrameCnt, Vid_FrameCnt)); 
  end
  % delete logger so that we don't accidentally keep recording into it.
  vid.DiskLogger = [];
end

%tdt: SET MODE TO RECORD
TDT_SetMode(xDA, 'standby'); % stop record if it's already running
[SysMode, recTank, recBlock, BlockStart_unixtime] = ...
  TDT_SetMode(xDA, 'record');

disp(sprintf('\nRecording started:\n  Tank: %s\n  Block: %s\n  T_start: %d\n',...
  recTank, recBlock, BlockStart_unixtime));

[recTankPath recTankName] = fileparts(recTank);
[topdir tanksdir] = fileparts(recTankPath);

%videoPath = [topdir filesep 'Video' filesep];
% Save video in TDT Data Tank folder for this block (makes copying data
% simpler)
videoPath = [recTank filesep recBlock filesep];

if ~exist(videoPath, 'dir'), mkdir(videoPath); end

%tdt: RESET FRAME COUNT
xDA.SetTargetVal([dev '.ResetFrameCnt'], 1);
xDA.SetTargetVal([dev '.ResetFrameCnt'], 0);

%vid: SET UP FILE LOGGER FOR VIDEO
vid.LoggingMode = 'disk';
diskLogger = VideoWriter([videoPath filesep recTankName '_' recBlock '_' datestr(now,30) '.mp4'], 'MPEG-4');
diskLogger.Quality = 95;
diskLogger.FrameRate = src.AcquisitionFrameRateAbs; % may be different than requested

disp(sprintf('Video file: %s\n',...
  [diskLogger.Path filesep diskLogger.Filename]));


%vid: clear old logger and replace with new one
delete(vid.DiskLogger);
vid.DiskLogger = diskLogger;

%vid: START VIDEO OBJECT
preview(vid);
start(vid);

% %tdt: TURN ON VIDACQ
xDA.SetTargetVal([dev '.VidAcq'], 1);

% % tdt: TURN ON VIDACQ if it was on when we started recording
% xDA.SetTargetVal([dev '.VidAcq'], VidOn);

disp(sprintf('   --> DONE\n'));
