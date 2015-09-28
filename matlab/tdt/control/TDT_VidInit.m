function [vid, src] = TDT_VidInit(format, fps, exposure_ms, binning, gain)
% Set up camera to synchronize with TDT (requires Image Acq toolbox)
%
% [vid src] = TDT_VidInit (format, fps, exposure_ms, binning, gain)
% 
% Sets up an AVT Manta G125C camera for behavior recording, including
% external trigger and digital outputs (one tick per frame) for TDT control
% of video recording, and synchronization of video with TDT hardware.
%
% Examples:
%  
% [vid src] = TDT_VidInit('BayerRG8',20,30,2,10);
% [vid src] = TDT_VidInit('Mono8',20,30,1,10);

%% defaults
if ~exist('fps', 'var') || isempty(fps), fps = 20; end
if ~exist('exposure_ms', 'var') || isempty(exposure_ms), exposure_ms = 40; end
if ~exist('binning', 'var') || isempty(binning), binning = 3; end
if ~exist('gain', 'var') || isempty(gain), gain = 10; end

% Set up videoinput object (uses Image Acquisition Toolbox), for use with
% AVT Manta G-125C monochrome camera

%% Reuse vid object if possible
vid = imaqfind('Tag', 'TDTVideo');
if iscell(vid),
  vid = [vid{:}];
end

switch numel(vid)
  case 0
    % create new videoinput
    vid = videoinput('gige', 1, format,'Tag','TDTVideo');
  case 1    
    %reuse videoinput
    if ~strcmp(vid.VideoFormat, format),
        delete(vid);
        vid = videoinput('gige', 1, format,'Tag','TDTVideo');
    end        
  otherwise
    error('More than 1 TDTVideo object already exists.')
end

%% Set up videoinput and source objects
vid.FramesPerTrigger = Inf;
vid.TriggerRepeat = 0;

src = getselectedsource(vid);
% make sure your network card is setup to support jumbo frames this big
src.PacketSize = 16000; 

%% Set up imaging parameters
src.AcquisitionFrameRateAbs = fps; %fps (may be modified by camera)
src.ExposureTimeAbs = exposure_ms * 1000; % microseconds. 40k = 1/50s
% src.BinningVertical = binning; % reduce frame size
% src.BinningHorizontal = binning;
src.DecimationVertical = binning; % reduce frame size
src.DecimationHorizontal = binning;
%vid.ROIPosition = [41 0 (vid.VideoResolution-[82 0])];
vid.ROIPosition = [0 0 (vid.VideoResolution)];
src.Gain = gain; % 0 - 32; 10 results in reasonably low noise


%% set up hardware triggering to respond to TDT digital output 
% (per Manta Camera Controls manual p.46)
% Camera is setup to go into free-running mode while Input 1 is high

triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');

set(src, 'AcquisitionStartTriggerMode','On')
set(src, 'AcquisitionStartTriggerSource','Line1')
set(src, 'AcquisitionStartTriggerActivation','LevelHigh')

set(src, 'AcquisitionEndTriggerMode','On')
set(src, 'AcquisitionEndTriggerSource','Line1')
set(src, 'AcquisitionEndTriggerActivation','LevelLow')

set(src, 'AcquisitionRecordTriggerMode','Off')

src.FrameStartTriggerMode = 'Off';
%src.FrameStartTriggerSource = 'Freerun';


%% set up camera frame sync output to be read by TDT

% Use 'FrameReadout' (time of frame transfer from CCD to camera memory)
% instead of 'Exposure' to avoid missed frames. (In 'Exposure' mode, when 
% frame rate is limited by exposure time, the 'low' period between frames 
% can be very short (~us), and can be ignored by TDT's digital input 
% debouncing.)
src.SyncOut1SyncOutSource = 'FrameReadout';
