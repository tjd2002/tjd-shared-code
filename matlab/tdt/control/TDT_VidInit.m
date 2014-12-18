function [vid, src] = TDT_VidInit(fps, exposure_ms, binning)

%% defaults
if ~exist('fps', 'var') || isempty(fps), fps = 20; end
if ~exist('exposure_ms', 'var') || isempty(exposure_ms), exposure_ms = 40; end
if ~exist('binning', 'var') || isempty(binning), binning = 3; end

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
    vid = videoinput('gige', 1, 'Mono8','Tag','TDTVideo');
  case 1    
    %reuse videoinput
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
src.BinningVertical = binning; % reduce frame size by 9x
src.BinningHorizontal = binning;

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

src.FrameStartTriggerSource = 'Freerun';


%% set up camera frame sync output to be read by TDT

% Use 'FrameReadout' (time of frame transfer from CCD to camera memory)
% instead of 'Exposure', because when frame rate is exposure time limited,
% the 'low' period between frames can be very short, and be missed by 
% TDT's digital input debouncing, causing frame timestamps to be missed.
src.SyncOut1SyncOutSource = 'FrameReadout';
