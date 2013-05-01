function [wfs] = NlxSpike2wf (varargin)
% Load several Neuralynx .nse files into a waveform (wf) struct
% Defaults to loading all *.nse files from the current directory
%
% Specify files to load using:
%  'NSEDir': folder containing NSE files (defaults to current dir)
%  'NSEFilenames': pattern to match for files to be loaded (default: '*.nse')
%                  if a cell array, files matching each pattern will be
%                  loaded
%     --OR--
%  'NSEFilePaths' full path to each file (cell array of strings)
%
%  'TimeWin': start/end of time to be imported (m x 2 array; default: all data)
%  'Name': name of cont struct (default: built from file paths)

% e.g. To load only Spike/SE1.nse and Spike/SE3.nse:
%
%  >> wf = NlxNSE2wf ('NSEDir', '/home/tjd/data/wallaceII_day3/Spike',...
%                            'NSEFilenames', {'SE1.nse' 'SE3.nse'})

% TODO:
% -handle file regexps rather than shell globbing?

%% parse and process inputs
p = inputParser;
p.addParamValue('NSEFilePaths', []);
p.addParamValue('NSEDir', '');
p.addParamValue('NSEFilenames', '*.nse');
p.addParamValue('TimeWin', [-Inf Inf], @isnumeric);
p.addParamValue('Name', [], @ischar);
p.parse(varargin{:});

% parse inputs
a = p.Results;

% If user specifies full paths to each file, use them
if ~isempty(a.NSEFilePaths),
    NSEFiles = a.NSEFilePaths;
    if ~iscell(NSEFiles),
        NSEFiles = {NSEFiles};
    end

    % come up with a fallback name for the result structure
    cdatName = cell2str(a.NSEFilePaths);
    
else % build up file names from specified directories and regexps
    
    if isempty(a.NSEDir),
        warning(['No ''NSEDir'' provided; using current working directory: ' pwd]); %#ok
        a.NSEDir = pwd;
    end
    
    if ~isdir(a.NSEDir)
        error(['Specified directory does not exist: ' a.NSEDir]);
    end

    % make sure there's a separator at the end
    a.NSEDir(end+1) = filesep;
    
    % convert string arg to cell array
    if ischar(a.NSEFilenames),
       a.NSEFilenames = {a.NSEFilenames};
    end

    NSEFiles = {};
    for j = 1:numel(a.NSEFilenames),
        dtmp = dir([a.NSEDir a.NSEFilenames{j}]);
        matches = {dtmp.name};
        for k = 1:numel(matches)
            NSEFiles = [NSEFiles; {fullpath([a.NSEDir matches{k}])}];
        end
    end
    
     
end % we should now have a cell array of full paths, 1 per file

if ~exist('NSEFiles', 'var') || isempty(NSEFiles),
  error('No filenames matched');
end

for j = 1:numel(NSEFiles),
    NSEFile = NSEFiles{j};
    
    % show progress
    disp(['Loading NSE file ' num2str(j) ' of ' num2str(numel(NSEFiles)) ...
        ' : ' NSEFile]);
    
    % load file into temporary struct
    wf(j) = Nlx2MatSpike(NSEFile,...

  
% 
%     NSE = NlxLoadNSE(NSEFile,...
%                      'TimeWin', a.TimeWin,...
%                      'DataUnits', 'mV', ...
%                      'TimeUnits', 'seconds', ...
%                      'UnwrapBuffers', true, ...
%                      'CorrectFilterDelay', true);

FieldSelection = [1 1 1 1 1]; % [timestamps Sc_Numbers Cell_numbers Params Data_points]
ExtractHeader = 1;
ExtractMode = 1; % 1-all, 2-range, 3-list, 4-timestamp range, 5-ts list
ModeArray = []; % blank for mode 1, 2 elements for range (mode 2 or 4);
                     % n elements for list (mode 3 or 5). indexes (modes 2 and
                     % 3) are zero-indexed

[se.spiketimes se.info.chan cs.info.fS_fromfile NumberValid cs.samples cs.info.rawheader] = ...
    Nlx2MatCSC(a.Filename, FieldSelection, ExtractHeader, ExtractMode, ModeArray);


%% detect some unusual problems

% error if file reports any invalid samples
if ~a.IgnoreInvalid
    if min(NumberValid) < size(cs.samples,1),
        error('Invalid data reported, use ''IgnoreInvalid'' to override this error');
    end
end

% error if reported recording channels change during file
if any(diff(cs.info.chan,[],2)),
    error('Reported recording channel changes within file.')
else
    % keep one copy of channel list
    cs.info.chan = cs.info.chan(:,1)';
end

% error if reported samplerate changes during file
if any(diff(cs.info.fS_fromfile,[],2)),
    error('Reported sampling rate changes within file.')
else
    % keep one copy of reported sampling rate
    cs.info.fS_fromfile = cs.info.fS_fromfile(1);
end

%% parse header and calculate some useful values

cs.info.header = NlxParseHeader(cs.info.rawheader);

% header ADBitVolts value has too few significant bits, recalculate it
cs.info.mVperADunit_fromheader = cs.info.header.ADBitVolts * 1000;
cs.info.mVperADunit_actual = cs.info.header.InputRange / ...
    (cs.info.header.ADMaxValue * 1000);

% get reported sampling rate
cs.info.fS_fromheader = cs.info.header.SamplingFrequency;


%% convert datatypes and scale as requested

% convert sample datatype as requested (to save memory, e.g.)
if ~strcmp(a.Datatype, class(cs.samples)), % only convert if necessary
    cs.samples = cast(cs.samples, a.Datatype);
end

% convert samples from digitizer units to Volts/mV as requested

switch a.DataUnits
    case 'mV'
        cs.samples = cs.samples .* cs.info.mVperADunit_actual;
    case 'V'
        cs.samples = cs.samples .* (cs.info.mVperADunit_actual / 1000);
end
cs.sampleunits = a.DataUnits;

% convert times from microseconds to seconds as requested
switch a.TimeUnits
    case 'seconds'
        cs.bufftimes = cs.bufftimes ./ 1e6;
    case 'microseconds'
        % no convert
    otherwise
        error('bad TimeUnits');            
end
cs.timeunits = a.TimeUnits;

% correct group delay introduced by digital filtering (see NOTES)
if a.CorrectFilterDelay
    filtdelay = ...
        (cs.info.header.DspLowCutNumTaps + ...
         cs.info.header.DspHighCutNumTaps - 1) ...
         / (2 .* CheetahfS);
    cs.bufftimes = cs.bufftimes - filtdelay;
end
  
end

                 
