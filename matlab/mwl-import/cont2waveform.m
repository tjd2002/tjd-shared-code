function wf = cont2waveform(cont, varargin)
% CONT2WAVEFORM - extracts waveform windows from a continuous eeg record. (can mimic Matt Wilson's recording software 'AD')
%
% cont2waveform(cont, ...
%  name/value args, any order
%
% one or both of:
%                'upthresh', [], 
%                'downthresh', [],
%
%                'rateperchan', [], 
%                'spikewin', 32,
%                'prespikepts', 6,
%                'postignore', 26,
%                'ignorefirstbuff', false,
%
% Defaults are designed to mimic Matt Wilson's AD acquisition software.
%
% cont is a struct (as returned by mwlIO load eeg file) with two members:
%  'timestamp' [1 x b uint32]
%  'data' [n x s x b int16]
%
%  where b = number of buffers
%       timestamps are for first sample of buffer
%       s = number of samples per buffer
%       n = number of channels (4 for tetrode, 8 for octrode, etc)  
%
% -upthresh and downthresh are thresholds in ad board units (i.e. eeg-file
% units, -4096 to 4096, int16). If both are specified, we trigger on
% upward and downward crossings of the two values. 'postignore' is enforced
% across the union of the crossings. (E.g. if you trigger on a downward
% crossing, you have to wait 'postignore' samples before you can trigger on
% either an upward or a downward crossing.)
%
% -rateperchan: sampling rate--only needed if only one buffer provided
%
% -spikewin is number of waveform points to extract per spike (AD = 32)
%
% -prespikepts is number of points before the threshold crossing to include
% (AD = 6)
%
% -postignore is the number of points after a threshold crossing where we
% ignore other threshold crossings (AD = 26).
%
% -'ignorefirstbuff' takes a logical value. To detect spikes that cross
% buffers, call cont2waveform with an extra 'previous' buffer and set this
% to true (default false)
%
% Tom Davidson (tjd@alum.mit.edu) 5/25/2006

% wishlist/minutiae
% -different thresholds per chan? (array of thresholds)
% -test for dropped buffers/breaks in eeg file (to avoid spurious threshold
%  xings across breaks?--just throw out spikes that cross 'bad' edges)

% struct cont has correct members, type, dims?
if ~(isstruct(cont) && ...
     all(ismember(fieldnames(cont), {'timestamp', 'data'})) && ...
     isnumeric(cont.timestamp) && isnumeric(cont.data) && ...
     ndims(cont.data)==3);
  error('Invalid data');
end

a = struct('upthresh',[], ...
           'downthresh',[],...
           'spikewin', 32, ...
           'prespikepts', 6, ...
           'postignore', 26, ...
           'rateperchan', [], ...
           'ignorefirstbuff',0);

a = parseArgs(varargin,a);

% validate input args:

bufflen = size(cont.data,2);
data = cont.data(:,:);

% calculate # timestamps per sample (use sampling rate, if provided.)
if ~isempty(a.rateperchan),

  % assume that timestamps are 0.1 msec
  tspersamp = 10000/rateperchan;

else
  if length(cont.timestamp) > 1
    tspersamp = double((cont.timestamp(end) - cont.timestamp(1))) / ...
                (size(data,2) - bufflen);
  else 
    error(['Can''t determine sampling rate from single buffer; provide as ' ...
           '''rateperchan'' argument.']);
  end
end

if isempty(a.upthresh) && isempty(a.downthresh),
  error('At least one of ''upthresh'' or ''downthresh'' must be provided.');
end

if (~isnumeric(a.upthresh) ||  ...
    ~isnumeric(a.downthresh)),
  error('Threshold values must be numeric.');
end

for argstring = {'spikewin' 'prespikepts' 'postignore'},
  argstring = argstring{:};
  if ~isnumeric(a.(argstring)),
    error(['Argument ''' argstring ''' must be numeric']);
  end
  if a.(argstring) < 0,
    error(['Argument ''' argstring ''' must be non-negative']);
  end
  if rem(a.(argstring),1),
    error(['Argument ''' argstring ...
           ''' must be an integer (# of samples)''']);
  end
end

if a.prespikepts + a.postignore > a.spikewin,
  error(['''prespikepts''(' num2str(a.prespikepts) ') +' ...
         ' ''postignore''(' num2str(a.postignore) ') is larger than ' ...
         'spike window size ''spikewin''(' num2str(a.spikewin) ')']);
end

postspikepts = a.spikewin - a.prespikepts;

% cont is at least a.spikewin samples
if size(data,2) < a.spikewin
	error('EEG data shorter than spike window');
end


% Look for threshold crossings: all channels below -> any channel above 
% (taking 'max' looks across all channels) (analogous operation using
% 'min' for downward going crossings).
%
% Note due to diff we find the first point *before* the crossing.
%
% Note we find 'kisses', when the signal is = threshold, but doesn't
% cross it.
%
% Ignore first 'prespikepts' and last 'postspikepts' points, so we don't
% read off the end of the array later. Note this means xings is indices 
% into an array spikewin smaller than data. Have to add prespikepts 
% to all indices when accessing, below

upxings = []; 
if ~isempty(a.upthresh),  
  mcth = max(data);
  % When diff of logical == 1, we have a xing from false to true.
  % If signal starts above thresh, not counted as xing.
  upxings = find(diff(mcth > a.upthresh) == 1);
end

% repeat for downward xings
downxings = [];
if ~isempty(a.downthresh),
  mcth = min(data);
  downxings = find(diff(mcth < a.downthresh) == 1);
end

% sort them both into one list
xings = sort([upxings downxings]);

% logical array of xings to ignore--assume none to start.
% (Faster than resizing xings for each exclusion)
xingsignore = false(length(xings),1);

% don't use points in last postspikepts (to avoid reading off end of
% buffer)
xingsignore(find(xings < size(data,2)-postspikepts, 1, 'last')+1:end) = true;

% don't use points in first prespikepts (to avoid reading off beginning of buffer)
xingsignore(1:find(xings > a.prespikepts, 1, 'first') - 1) = true;
	
%%% postignore
%
% choose only crossings that are at least 'postignore' apart to,
% e.g. avoid creating a new spike for each channel if they cross at
% different times.
%
% NOTE: we can't just throw away all points after short gaps, since the
% postignore gap is only 'renewed' by a threshold xing that results in a
% spike being grabbed, and *not* by an ignored spike.

if length(xings > 1),

  % nextoksamp is the next sample after the ignore window. 
  % init to 0 so 1st spike will be good.
  nextoksamp = 0; 
  
  % iterate over xings
  for j = 1:length(xings),
    
    if xings(j) > nextoksamp,
      % keep this one as ok, ignore xings in next postignore samps
      nextoksamp = xings(j) + a.postignore;
      
    else 
      % ignore this xing, don't update nextoksamp
      xingsignore(j) = true;
      
    end
  end
end


%%% ignorefirstbuffer
%
% handle ignoring first buffer (so that we catch spikes on the edge of
% buffer). Do this after calc'ing spikes for first buffer to correctly
% handle postignore for spikes in the 'previous' buffer.

if a.ignorefirstbuff;
  xingsignore(1:find(xings >= bufflen-postspikepts, 1, 'first') - 1) = true;
end


%%% throw out ignored xings:

xings = xings(~xingsignore);
  

% initialize output struct
nspikes = size(xings,2);
wf = struct('timestamp',[],'waveform',[]);
wf.timestamp = zeros(1,nspikes,class(cont.timestamp));
wf.waveform = zeros(size(data,1),a.spikewin,nspikes, class(data));

% fill output struct
for j = 1:nspikes ;
  
	% lookup timestamp for this buffer (keep sync'd with AD's clock)
	% (+1 is because timestamp array is 1-indexed)
	wf.timestamp(j) = cont.timestamp(floor(xings(j)/bufflen)+1) + rem(xings(j),bufflen)*tspersamp;

	% grab spikewin pts around threshold xing
	wf.waveform(:,:,j) = data(:,xings(j)-a.prespikepts:xings(j)+postspikepts-1);
end
