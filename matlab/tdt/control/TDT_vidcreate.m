%%% DEPRECATED--use TDT_vid_expt_script

%% TODO
% -parameterize video/tank storage locations
% -write functions/dedupe code
%  -ToStandby/Idle (close video files if exist)
%  -ToPreview
%  -ToRecord
% -Confirm that dev.FrameCnt == #frames in video file
% -wrapper function for SetTargetVal 
%   -confirm device is running/accepting input
%   -check return values
%   -optionally 'toggle' a logical
%   -check target type/size?
% -Block to go to Idle instead of Preview
% -Go to 'Idle' instead of 'Standby'? (Will auto-reset frame count)
% -Configure/check whether variables are Stored/Fetched/Disabled (in idle
%  mode)
% -When stopping video record, explicitly check for arrival of last frame
% (instead of just waiting 0.5s)

%% MAYBE
% -Get data from tank?

%% DONE
% -add tank name/block name to video file name
% -Display tank/file names 

%% DEPENDENCIES/ASSUMPTIONS
% 1 Device connected (index 0), with paramater tags:
%  -VidAcq: Logical to start/stop video acquisition, controls digital out to
%   free-running video camera
%  -ResetFrameCnt: Logical to reset FrameCnt to 1 (toggle 1/0 to reset)

%% INITIALIZE TDT AND CAMERA
[xDA] = TDT_Init;

%SETUP TANK LOCATION (ADDTANK?)

%[vid, src] = TDT_VidInit('Mono8', 20,40,3,0); % format, fps, exposure, binning
vid = TDT_VidInit('BayerRG8', 20,20,1,10); % format, fps, exposure, binning, gain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ->Preview (no logging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n\n*** GOING TO PREVIEW MODE ***'));

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
end

%tdt: SET MODE TO PREVIEW
TDT_SetMode(xDA, 'standby'); % stop preview if it's already running
TDT_SetMode(xDA, 'preview');

%tdt: RESET FRAME COUNT
xDA.SetTargetVal([dev '.ResetFrameCnt'], 1);
xDA.SetTargetVal([dev '.ResetFrameCnt'], 0);

%vid: START VIDEO PREVIEW
preview(vid);

%tdt: TURN ON VIDACQ
xDA.SetTargetVal([dev '.VidAcq'], 1);
xDA.GetTargetVal([dev '.VidAcq'])

disp(sprintf('   --> DONE\n'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ->Record into new block, new video file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n\n*** GOING TO RECORD MODE ***'));

%tdt: Save VIDACQ state
VidOn = xDA.GetTargetVal([dev '.VidAcq']);

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
videoPath = [topdir filesep 'Video' filesep];
mkdir(videoPath);

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
% xDA.SetTargetVal([dev '.VidAcq'], 1);

% tdt: TURN ON VIDACQ if it was on when we started recording
xDA.SetTargetVal([dev '.VidAcq'], VidOn);

disp(sprintf('   --> DONE\n'));
