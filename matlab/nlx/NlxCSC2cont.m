function [cdat_lfps] = NlxCSC2cont (varargin)
% Load all CSC files from a recording session into a single cont struct
% Defaults to loading all *.ncs files from the directory called 'LFP'
%
% Specify files to load using:
%  'CSCFilePaths' full path to each CSC file (cell array of strings)
%   -OR-
%  'SessDir' : top-level folder (string; defaults to current dir)
%  'CSCDir' : folder containing CSC files (defaults to 'LFP')
%  'CSCFilenames' : pattern to match for files to be loaded (default: '*.ncs')
%               if a cell array, files matching each pattern will be
%              loaded
%
% e.g. To load only LFP/LFP1.ncs and LFP/LFP3.ncs:
%
%  >> cdat_LFP = NlxCSC2cont('SessDir', 'wallaceII_day3',...
%                               'CSCDir', 'LFP',...
%                               'CSCRegexp', {'LFP1*' 'LFP3*'})

% TODO:
% -handle file regexps

%% parse and process inputs
p = inputParser;
p.addParamValue('CSCFilePaths', []);
p.addParamValue('SessDir', [], @isdir);
p.addParamValue('CSCDir', 'LFP');
p.addParamValue('CSCFilenames', '*.ncs');
p.addParamValue('TimeWin', [-Inf Inf], @isnumeric);
p.addParamValue('Name', [], @ischar);
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
    
    if isempty(a.SessDir),
        warning(['using current working directory for SessDir: ' pwd]);
        a.SessDir = pwd;
    end
    
    % get our working directory
    wd = [a.SessDir '/' a.CSCDir];
    if ~isdir(wd)
        error(['Specified directory does not exist: ' wd]);
    end
    
    if ischar(a.CSCFilenames),
       a.CSCFilenames = {a.CSCFilenames};
    end

    CSCfiles = {};
    CSCdir = [a.SessDir '/' a.CSCDir '/'];
    for j = 1:numel(a.CSCFilenames),
        dtmp = dir([CSCdir a.CSCFilenames{j}]);
        matches = {dtmp.name};
        for k = 1:numel(matches)
            CSCfiles = [CSCfiles; {fullpath([CSCdir matches{k}])}];
        end
    end
    
    % come up with a fallback name for the result structure
    if numel(a.CSCFilenames) == 1
        CSCFilenamesStr = a.CSCFilenames{1};
    else
        CSCFilenamesStr = cell2str(a.CSCFilenames);
    end
    cdatName = [CSCdir CSCFilenamesStr];
     
end % we should now have a cell array of full paths, 1 per file


for j = 1:numel(CSCfiles),
    CSCfile = CSCfiles{j};
    
    % show progress
    disp(['Loading CSC file ' num2str(j) ' of ' num2str(numel(CSCfiles)) ...
        ' : ' CSCfile]);
    
    % load file into temporary struct
    csc = NlxLoadCSC(CSCfile,...
                     'TimeWin', a.TimeWin,...
                     'DataUnits', 'mV', ...
                     'TimeUnits', 'seconds', ...
                     'UnwrapBuffers', true, ...
                     'CorrectFilterDelay', true);
    
    % convert to contstruct (compatible with Tom's viewer and cont* functions)
    cdat_lfps{j} = imcont('neuralynxCSC', csc);

end

% merge multiple channels into one contstruct
if numel(cdat_lfps)>1
    cdat_lfps = contcombine(cdat_lfps{1}, cdat_lfps(2:end));
end

% assign name to struct
if ~isempty(a.Name)
    cdat_lfps.name = a.Name;
else
    warning (['No ''Name'' provided; using default struct name: ' cdatName]);
    cdat_lfps.name = cdatName;
end
                 