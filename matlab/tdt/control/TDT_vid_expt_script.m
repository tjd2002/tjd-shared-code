%% TODO
% -parameterize video/tank storage locations
% -wrapper function for SetTargetVal 
%   -confirm device is running/accepting input
%   -check return values
%   -optionally 'toggle' a logical
%   -check target type/size?
% -Watch (and optionally set) experimental parameters every 1sec: FP frequencies, LED
% power, filters, LFP ref channel, etc. Save metadata in a .mat struct
% alongside each recorded block, with a new struct each time a setting
% changes?

%% MAYBE
% -Do FP filtering in Matlab?
% -Have Matlab watch the VidAcq state, and open video recording whenever 
%  it changes.(i.e. make TDT the primary controller).
% -Get data from tank?
% -Configure/check whether variables are Stored/Fetched/Disabled (in idle
%  mode)
% -When stopping video record, explicitly check for arrival of last frame
% (instead of just waiting 0.5s)

%% DONE
% -add tank name/block name to video file name
% -Display tank/file names 
% -write functions/dedupe code
%  -ToStandby/Idle (close video files if exist)
%  -ToPreview
%  -ToRecord
% -Confirm that dev.FrameCnt == #frames in video file
% -Block to go to Idle instead of Preview

%% DEPENDENCIES/ASSUMPTIONS
% 1 Device connected (index 0), with paramater tags:
%  -VidAcq: Logical to start/stop video acquisition, controls digital out to
%   free-running video camera
%  -ResetFrameCnt: Logical to reset FrameCnt to 1 (toggle 1/0 to reset)

%% INITIALIZE TDT AND CAMERA
[xDA] = TDT_Init;
%SETUP TANK LOCATION (ADDTANK?)
%[vid, src] = TDT_VidInit('Mono8', 20,40,3,0); % format, fps, exposure, binning
[vid src] = TDT_VidInit('BayerRG8', 20,50,2,20); % format, fps, exposure (ms), binning, gain

%% Preview / Stop recording
TDTvid_ToPreview(xDA,vid);

%% Record
TDTvid_ToRecord(xDA,vid);

%% Idle
TDTvid_ToIdle(xDA,vid);

