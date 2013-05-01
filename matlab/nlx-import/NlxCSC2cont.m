function [cdat_cscs] = NlxCSC2cont (varargin)
% Load several Neuralynx CSC files into a single cont struct
% Defaults to loading all *.ncs files from the current directory
%
% Specify files to load using:
%  'CSCDir': folder containing CSC files (defaults to current dir)
%  'CSCFilenames': pattern to match for files to be loaded (default: '*.ncs')
%                  if a cell array, files matching each pattern will be
%                  loaded
%     --OR--
%  'CSCFilePaths' full path to each CSC file (cell array of strings)
%
%  'TimeWin': start/end of time to be imported (m x 2 array; default: all data)
%  'Name': name of cont struct (default: built from file paths)

% e.g. To load only LFP/LFP1.ncs and LFP/LFP3.ncs:
%
%  >> cdat_LFPs = NlxCSC2cont('CSCDir', '/home/tjd/data/wallaceII_day3/LFP',...
%                            'CSCFilenames', {'LFP1.ncs' 'LFP3.ncs'})

% TODO:
% -handle file regexps rather than shell globbing?

%% parse and process inputs
p = inputParser;
p.addParamValue('CSCFilePaths', []);
p.addParamValue('CSCDir', '');
p.addParamValue('CSCFilenames', '*.ncs');
p.addParamValue('TimeWin', [-Inf Inf], @isnumeric);
p.addParamValue('Name', [], @ischar);
p.addParamValue('Force_ignore_ts_errors', false, @islogical);
p.addParamValue('ts_syn_linmode', 'regress', @ischar);
p.parse(varargin{:});

% parse inputs
a = p.Results;

% If user specifies full paths to each file, use them
if ~isempty(a.CSCFilePaths),
    CSCFiles = a.CSCFilePaths;
    if ~iscell(CSCFiles),
        CSCFiles = {CSCFiles};
    end

    % come up with a fallback name for the result structure
    cdatName = cell2str(a.CSCFilePaths);
    
else % build up file names from specified directories and regexps
    
    if isempty(a.CSCDir),
        warning(['No ''CSCDir'' provided; using current working directory: ' pwd]); %#ok
        a.CSCDir = pwd;
    end
    
    if ~isdir(a.CSCDir)
        error(['Specified directory does not exist: ' a.CSCDir]);
    end

    % make sure there's a separator at the end
    a.CSCDir(end+1) = filesep;
    
    % convert string arg to cell array
    if ischar(a.CSCFilenames),
       a.CSCFilenames = {a.CSCFilenames};
    end

    CSCFiles = {};
    for j = 1:numel(a.CSCFilenames),
        dtmp = dir([a.CSCDir a.CSCFilenames{j}]);
        matches = {dtmp.name};
        for k = 1:numel(matches)
            CSCFiles = [CSCFiles; {fullpath([a.CSCDir matches{k}])}];
        end
    end
    
    % come up with a fallback name for the result structure
    if numel(a.CSCFilenames) == 1
        CSCFilenamesStr = a.CSCFilenames{1};
    else
        CSCFilenamesStr = cell2str(a.CSCFilenames);
    end
    cdatName = [a.CSCDir CSCFilenamesStr];
     
end % we should now have a cell array of full paths, 1 per file

if ~exist('CSCFiles', 'var') || isempty(CSCFiles),
  error('No filenames matched');
end

for j = 1:numel(CSCFiles),
    CSCFile = CSCFiles{j};
    
    % show progress
    disp(['Loading CSC file ' num2str(j) ' of ' num2str(numel(CSCFiles)) ...
        ' : ' CSCFile]);
    
    % load file into temporary struct
    csc = NlxLoadCSC(CSCFile,...
                     'TimeWin', a.TimeWin,...
                     'DataUnits', 'mV', ...
                     'TimeUnits', 'seconds', ...
                     'UnwrapBuffers', true, ...
                     'CorrectFilterDelay', true);
    
    if a.Force_ignore_ts_errors,
      warning('!!! Using permissive timestamp options, check max_tserr after loading!');
    end
    
    % convert to contstruct (compatible with Tom's viewer and cont* functions)    
    cdat_cscs{j} = imcont('neuralynxCSC', csc, 'ts_syn_linmode', a.ts_syn_linmode, ...
                          'ts_permissive', a.Force_ignore_ts_errors);

end

% merge multiple channels into one contstruct
if numel(cdat_cscs)>1
    cdat_cscs = contcombine(cdat_cscs{1}, cdat_cscs(2:end));
else
  cdat_cscs = cdat_cscs{1};
end

% assign name to struct
if ~isempty(a.Name)
    cdat_cscs.name = a.Name;
else
    warning (['No ''Name'' provided; using default struct name: ' cdatName]);
    cdat_cscs.name = cdatName;
end
                 
