%% tdt2mat.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import data from TDT system recording into Matlab structure
%
% S = tdt2mat(filepath, tank, blk, event)
%
% INPUTS
%   filepath: folder where tank is stored
%   tank: tank name
%   blk: block name
%   event: event name (ie 'Wave')
%
% OUTPUT
% data structure with fields
%     sampling_rate: (Hz)
%        timestamps: timestamp corresponding to each data row
%          channels: channel corresponding to each data row
%              data: data in rows, npoints columns
%           npoints: number of data points per timestampt
%            nchans: number of channels
%
%
% SOURCE: http://jaewon.mine.nu/jaewon/2010/10/04/how-to-import-tdt-tank-into-matlab/
%
%
% Example usage:
%{
filepath = 'C:\TDT\OpenEx\Tanks';
tank = 'SETUPTANK';
blk = 'Test1-2';
event = 'Wave';

S = tdt2mat2(filepath, tank, blk, event)
%}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TODO:
% -Handle epoc_store secondary stores (what format code? where is the
% non-strobe data stored???)
% -Handle SlowStores (same issue as epoc_store?) OpenEx DataList internally
% -Select channels/timewins for loading (avoid out-of-memory)


function S = tdt2mat(filepath, tank, blk, event)



filename = strcat(tank, '_', blk);

tev_path = [filepath filesep tank filesep blk filesep filename '.tev'];
tsq_path = [filepath filesep tank filesep blk filesep filename '.tsq'];




% open the files
tev = fopen(tev_path);

% count number of tsq records (40 bytes/record)
tsq = fopen(tsq_path); fseek(tsq, 0, 'eof'); ntsq = ftell(tsq)/40; fseek(tsq, 0, 'bof');

% read from tsq
data.size      = fread(tsq, [ntsq 1], 'int32',  36); fseek(tsq,  4, 'bof');
data.type      = fread(tsq, [ntsq 1], 'int32',  36); fseek(tsq,  8, 'bof');
data.name(:,1) = fread(tsq, [ntsq 1], '*char*1', 39); fseek(tsq, 9, 'bof');
data.name(:,2) = fread(tsq, [ntsq 1], '*char*1', 39); fseek(tsq, 10, 'bof');
data.name(:,3) = fread(tsq, [ntsq 1], '*char*1', 39); fseek(tsq, 11, 'bof');
data.name(:,4) = fread(tsq, [ntsq 1], '*char*1', 39); fseek(tsq, 12, 'bof');
data.chan      = fread(tsq, [ntsq 1], 'ushort', 38); fseek(tsq, 14, 'bof');
data.sortcode  = fread(tsq, [ntsq 1], 'ushort', 38); fseek(tsq, 16, 'bof');
data.timestamp = fread(tsq, [ntsq 1], 'double', 32); fseek(tsq, 24, 'bof');
data.fp_loc    = fread(tsq, [ntsq 1], 'int64',  32); fseek(tsq, 24, 'bof');
data.strobe    = fread(tsq, [ntsq 1], 'double', 32); fseek(tsq, 32, 'bof');
data.format    = fread(tsq, [ntsq 1], 'int32',  36); fseek(tsq, 36, 'bof');
data.frequency = fread(tsq, [ntsq 1], 'float',  36);

% TDT timestamps are in Unix time (seconds since 1/1/1970--
% http://en.wikipedia.org/wiki/Unix_time )


allnames = unique(data.name,'rows');
%exclude 'names' that contain ASCII Null character, these are special TDT
%codes.
allnames = allnames(all(double(allnames)~=0,2),:);


% Read A/D samples from tev
row = ismember(data.name, event, 'rows');

if sum(row) == 0
  disp(['Requested store name ''' event ''' not found (case-sensitive).'])
  disp('File contains the following TDT store names:');
  disp(allnames);
  error('TDT store name not found.');
end

table = { 'float',  1, 'float';
          'long',   1, 'int32';
          'short',  2, 'short';
          'byte',   4, 'schar'; }; % a look-up table
first_row = find(row,1);
S.format    = data.format(first_row)+1; % from 0-based to 1-based

S.storename = event;
S.sampling_rate = data.frequency(first_row);
S.timestamps    = data.timestamp(row);
S.channels      = data.chan(row);


fp_loc  = data.fp_loc(row);

if S.format ~=5
  nsample = (data.size(first_row)-10) * table{S.format,2};
  S.data = zeros(length(fp_loc),nsample);
  for n=1:length(fp_loc)
    fseek(tev,fp_loc(n),'bof');
    S.data(n,1:nsample) = fread(tev,[1 nsample],table{S.format,3});
  end
  S.npoints = nsample;
else % epoc_stores and slow_stores have data as a float in the 'strobe' field.
  S.data = data.strobe(row);
  S.npoints = 1;
  % epoc_store events list all strobe data as channel 0; 
  % slow_store events use channel numbering
  S.channels = data.chan(row);
  if all(S.channels == 0),
      S.channels = repmat(1, size(S.channels));
  end
end

% Which channels are present in the file? 
S.chans = sort(unique(S.channels))';

