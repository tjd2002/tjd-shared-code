function S = TDT_Import(filepath, tank, blk, name, chans, timewin,...
    return_timestamps);
%% TDT_Import.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import data from TDT system recording into Matlab structure
%
% S = TDT_Import(filepath, tank, blk, event, chans)
%
% INPUTS
%   filepath: folder where tank is stored
%    tank: tank name
%     blk: block name
%    name: store/event name (ie 'Wave')
%
% OPTIONAL INPUTS
%   chans: for 'Wave' events, array of channel #s to import 
%          (Optional. []/default = import all chans)
% timewin: 1x2 array of start/end time to import (NOT CURRENTLY SUPPORTED)
% return_timestamps: 
%          flag whether to populate the 'timestamps' field with the time of
%          each sample. Redundant with tstart/tend/samplerate, so defaults
%          to false.
%
% OUTPUT
% data structure with fields
%         storename: 4-character Store name
%              data: data in rows, 1 column per channel
%             chans: TDT channel numbers, 1 per data column
%            tstart: start time of extracted data, relative to when recording started (s) 
%              tend: end time of extracted data, relative to when recording started (s) 
%     sampling_rate: (Hz)
%       format_code:
%
% derived from: 
%  http://jaewon.mine.nu/jaewon/2010/10/04/how-to-import-tdt-tank-into-matlab/
%
% Example usage:
%{
filepath = 'C:\TDT\OpenEx\Tanks';
tank = 'SETUPTANK';
blk = 'Test1-2';
event = 'Wave';
chans = [2 4]; % extract only 2 channels

S = TDT_Import(filepath, tank, blk, event, chans)
%}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TDT Import notes
% -'Tick' is not every second. Seems to be issued every n samples, where n
%  is  round(sampling_rate)-1. Could be floor instead of round. For 12k
%  (12207.03125 Hz) system clock, 'Tick' is issued every 12206 samples =
%  0.99991552 s.) Take-home: don't rely on 'Tick's for timing.
% -Sampling clock is same as timestamp clock, so don't have to worry about
%  skew, drift, nonlinearities. Woohoo! 
%
% TODO:
% -How to deal with outputting different data types? TDT uses substructures 
%  for each data type (Wave, Strobe, etc)
% -Select channels/timewins *during loading* (to avoid out-of-memory conditions)
% -Don't load all data at first
% -Support other TDT data types: snippets, etc.
% -Check out TDT-provided SEV2mat.m (which reads .sev files, and doesn't use ActiveX!)
% -'info' field for true start time, etc?
% -Handle intermittently recorded Waves (using timestamps field)
% -Handle epoc_store secondary stores (what format code? where is the
% non-strobe data stored???)
% -Handle SlowStores (same issue as epoc_store?) OpenEx DataList internally
% -convert formats to string. (OpenDeveloper manual EvTypeToString)
% -separate function to get list of store names
%
% DONE:
% -save a lot of memory by having fread use correct Matlab types (single,
%  uint32, etc.) instead of double for everything. Both .tsq and .tev
% -un-bufferize wave data
% -test for timestamp errors in buffers (no need to interpolate timestamps
%  on each buffer first
% -rename to avoid confusion with TDT2mat: TDT_LoadStore, TDT_ParseHeader?


% handle missing/optional inputs
if ~exist('name', 'var'), name= []; end; % will result in list of available store names being displayed
if ~exist('chans', 'var'), chans = []; end;
if ~exist('timewin', 'var'), timewin= []; end;
if ~exist('return_timestamps', 'var'), return_timestamps = false; end;

% Haven't implemented 'timewin' yet
if ~isempty(timewin),
    error('''timewin'' argument not currently supported; omit or leave empty to import all data.');
end

filename = strcat(tank, '_', blk);

tev_path = [filepath filesep tank filesep blk filesep filename '.tev'];
tsq_path = [filepath filesep tank filesep blk filesep filename '.tsq'];

S = struct(...
    'storename', [],...
    'data', [],...
    'chans', [],...
    'tstart', [],...
    'tend', [],...
    'timestamps', [],...
    'sampling_rate', [],...
    ...
    'format_code', [],...
    'buff_timestamps', [],...
    'buff_channel', [],...
    'buff_data', [],...
    'buff_npoints', [],...
    't_rec_start_UnixTime', [],...
    't_rec_start_UTC', [],...
    'max_ts_err', []);

% open the files
tev = fopen(tev_path);
tsq = fopen(tsq_path); 

recsize_bytes = 40;
% count number of tsq records (40 bytes/record)
fseek(tsq, 0, 'eof'); ntsq = ftell(tsq)/recsize_bytes; 

% read from tsq
fseek(tsq, 0, 'bof');
data.size      = fread(tsq, [ntsq 1], '*int32',  recsize_bytes-4); 

fseek(tsq, 4, 'bof');
data.type      = fread(tsq, [ntsq 1], '*int32',  recsize_bytes-4); 

fseek(tsq, 8, 'bof');
data.name      = fread(tsq, [4 ntsq], '4*uchar=>char', recsize_bytes-4)'; %reads column-wise

fseek(tsq, 12, 'bof');
data.chan      = fread(tsq, [ntsq 1], '*ushort', recsize_bytes-2);

% % 'sortcode' not currently used for import
% fseek(tsq, 14, 'bof'); 
% data.sortcode  = fread(tsq, [ntsq 1], '*ushort', recsize_bytes-2); 

fseek(tsq, 16, 'bof'); 
data.timestamp = fread(tsq, [ntsq 1], '*double', recsize_bytes-8); 

% this position (24) is either a pointer into .tev, or a double float as data 

fseek(tsq, 24, 'bof');
data.fp_loc    = fread(tsq, [ntsq 1], '*int64',  recsize_bytes-8);

% read 'strobe' below only if requested
fseek(tsq, 24, 'bof');
data.strobe    = fread(tsq, [ntsq 1], '*double', recsize_bytes-8); 

fseek(tsq, 32, 'bof');
data.format    = fread(tsq, [ntsq 1], '*int32',  recsize_bytes-4); 

fseek(tsq, 36, 'bof');
data.frequency = fread(tsq, [ntsq 1], '*float',  recsize_bytes-4);

allnames = unique(data.name,'rows');
% exclude 'names' that contain ASCII Null character (==0), these are 
% special TDT codes, not real stores.
allnames = allnames(all(double(allnames)~=0,2),:);

% Get rows matching requested store name
if isempty(name), 
    row = [];
else
    row = ismember(data.name, name, 'rows');
end

if sum(row) == 0
  disp(['Requested store name ''' name ''' not found (case-sensitive).'])
  disp('File contains the following TDT store names:');
  disp(allnames);
  error('TDT store name not found.');
end

% Second entry in the .tsq file always contains the timestamp of the start 
% of the block. TDT timestamps are in Unix time (seconds elapsed since 
% January 1, 1970 -- http://en.wikipedia.org/wiki/Unix_time )
S.t_rec_start_UnixTime = data.timestamp(2);
S.t_rec_start_UTC = [datestr(S.t_rec_start_UnixTime/86400 + datenum(1970,1,1)) 'Z'];

% See OpenDeveloper manual 'GetNPer' and 'DFromToString' (a typo for Data FORMat To String)
% Format codes are 0-5 in 6 rows below
%       { TDT type,  bytes, fread code,       Matlab type }
table = { 'float',   4,     'float=>single',  'single'; % 0 
          'long',    4,     'int32=>int32',   'int32';  % 1
          'short',   2,     'short=>int16',   'int16';  % 2
          'byte',    1,     'schar=>int8',    'int8';   % 3
          'double',  8,     'double=>double', 'double'; % 4 
          'qword'    8,     'int64=>int64',   'int64'}; % 5 not seen in the wild

first_row = find(row,1);
S.format_code    = data.format(first_row); 
format_idx = S.format_code+1; % from 0-based to 1-based for table lookup

S.storename = name;
S.sampling_rate = double(data.frequency(first_row));
S.buff_timestamps    = data.timestamp(row);
S.buff_channel      = data.chan(row);


fp_loc  = data.fp_loc(row);

if S.format_code ~=4
  nsample = ((4 .* data.size(first_row)) - 40) ./ table{format_idx,2};
  if rem(nsample,1), error('Non-integer number of samples--problem with size calculation'); end
  S.buff_data = zeros(length(fp_loc),nsample, table{format_idx,4});
  for n=1:length(fp_loc)
    fseek(tev,fp_loc(n),'bof');
    S.buff_data(n,1:nsample) = fread(tev,[1 nsample],table{format_idx,3});
  end
  S.buff_npoints = double(nsample);
else % epoc_stores and slow_stores have data as a float in the 'strobe' field.
  S.buff_data = data.strobe(row);
  S.buff_npoints = 1;
  % epoc_store events list all strobe data as channel 0; 
  % slow_store events use channel numbering
  S.buff_channel = data.chan(row);
  if all(S.buff_channel == 0),
      S.buff_channel = repmat(1, size(S.buff_channel));
  end
end

% Which channels are present in the file? 
S.chans = sort(unique(S.buff_channel))';

% 1: continuously-sampled 'Wave' store
% Can we unbufferize other formats?
if S.format_code == 0

    % validate requested channels
    if isempty(chans) % default to all channels, sorted
      chans = S.chans; 
    elseif ~all(ismember(chans, S.chans)),
      % error: display bad values
      S.chans
      chans
      error('Requested ''chans'' channels not present in file');
    end
    
    nchans = numel(chans);

    % how many buffers per channel? (calculate for first channel)
    chanione = S.buff_channel==chans(1);
    nbuffs = sum(chanione);
    
    % initialize 3-D array (bufferno, bufferlen, chan)
    dat = zeros(nbuffs, S.buff_npoints, nchans, class(S.buff_data));
    
    % iterate over all channels
    for k = 1:nchans,
        % keep user order of channels requested
        chan = chans(k);
        chani = S.buff_channel==chan;
        if sum(chani) ~= nbuffs,
            error('Mismatch in number of buffers for channel %d', chan)
        end
        
        % get (still buffered) data for each channel
        dat(:,:,k) = S.buff_data(chani,:);
    end
    
    % reshape buffered data into single columns, one per channel
%     datp = permute(dat,[2 1 3]);
%     datpr = reshape(datp,[],nchans);
%     dat = datpr;    
    S.data = reshape(permute(dat,[2 1 3]), [], nchans);
    
    % Process timestamps
    ts = S.buff_timestamps(chanione);
    % convert from TDT timestamps to 'seconds from block start'
    ts = ts-S.t_rec_start_UnixTime; 
    
    % Check for missing timestamps (one per buffer)
    ts_syn = linspace(ts(1),ts(end), numel(ts))'; % 'synthetic' timestamps
    S.max_ts_err = max(abs(ts_syn-ts));
    
    if S.max_ts_err > 1./S.sampling_rate;
        error('Timestamp error: at least one sample gap or repeat');
    end

    S.tstart = ts(1);
    S.tend = ts(end)+(S.buff_npoints-1)/S.sampling_rate;
    
    % Interpolate timestamps over entire file (more accurate than doing it
    % buffer-wise, since sample clock is exact).
    if return_timestamps
        S.timestamps = ...
            linspace(ts(1),...
            ts(end)+(S.buff_npoints-1)/S.sampling_rate,...
            numel(S.data));
    end
    
    % discard some fields we don't want to return for 'wave's
    S = rmfield(S, {'buff_npoints', 'buff_channel', 'buff_data', 'buff_timestamps'});
end
