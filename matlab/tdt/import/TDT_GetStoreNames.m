function [names types] = TDT_GetStoreNames(filepath, tank, block)
% TDT_GetStoreNames - Get names and types of all 'stores' in a TDT recording
%
% [names types] = TDT_GetStoreNames(filepath, tank, blk)
%
% INPUTS
%   filepath: folder where tank is stored
%       tank: tank name
%      block: block name
%
% OUTPUT
%   names: array of 4-character store names
%
% Example usage:
%
% filepath = 'C:\TDT\OpenEx\Tanks';
% tank = 'SETUPTANK';
% blk = 'Test1-2';
%
% names = TDT_GetStoreNames(filepath, tank, blk);
%
% Author: Tom Davidson - tjd@alum.mit.edu - Dec. 2014
% Based on TDT format info from: http://jaewon.mine.nu/jaewon/2010/10/04/how-to-import-tdt-tank-into-matlab/

% TODO:

%% Constants and data types

% .tsq file record size, in bytes
tsq_recsize = 40;

% TDT codes for event types:
TDT_event_table = ...
   ...code   name    is_a_store?   comment
   {    0 'UNKNOWN'   0;...      % first event in TSQ file (0x00000000)
    34817 'MARK'      0;...      % start/end of recording (0x00008801)
    32768 'HAS_DATA'  0;...      % not sure when used (0x00008000)
      257 'Strobe+'   1;...      % strobe event onset (0x00000101)
      258 'Strobe-'   1;...      % strobe event offset (rarely used; 0x00000102)
      513 'Scalar'    1;...      % scalar list (unused?; 0x00000201)
    33025 'Stream'    1;...      % continuous data (0x00008101)
    33281 'Snip'      1;...      % 'snippets' of continuous data (0x00008201)
    };        
%   33073 'Stream'       1       % Mentioned in TDT's TDT2mat.m as a code
                                 % for RS4 headers (would be 0x0008131), 
                                 % detected using bitand(type, 33025 (0x8101)
                                 % and set to 'Stream' type.

% Get codes of TDT-generated events, with no corresponding store name.
rowidx = ~vertcat(TDT_event_table{:,3});
TDT_nonstore_codes = [TDT_event_table{rowidx,1}];

%% Open .tsq file containing data headers
filename = strcat(tank, '_', block);
tsq_path = [filepath filesep tank filesep block filesep filename '.tsq'];
tsq = fopen(tsq_path);
if tsq == -1
    error('Could not open file ''%s'' for reading.', tsq_path);
end

%% get number of records (40 bytes/record)
% first 'size' record contains # of bytes expected in file
frewind(tsq);
ntsq = fread(tsq, 1, '*int32') ./ tsq_recsize;

%% read Store name and type for all data records in file
fseek(tsq, 8, 'bof'); % 8-byte offset for 'name' field
data.name      = fread(tsq, [4 ntsq], '4*uchar=>char', tsq_recsize-4)'; %reads column-wise

fseek(tsq, 4, 'bof');
data.evtype      = fread(tsq, [ntsq 1], '*int32',  tsq_recsize-4); 

% exclude headers for special TDT-generated events (not real stores)
validi = true(size(data.evtype));
for nonstore_code = TDT_nonstore_codes,
    validi = validi & data.evtype ~= nonstore_code;
end

% return list of store names occurring at least once in block
[names idx] = unique(data.name(validi,:), 'rows');

evtype_valid = data.evtype(validi);
evtypes = evtype_valid(idx);
%types = TDT_EvTypeToString(evtypes);

