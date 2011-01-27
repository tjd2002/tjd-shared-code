function wf = cdat2waveform(c, varargin)
% CDAT2WAVEFORM - extracts waveform windows from a continuous eeg record. (can mimic Matt Wilson's recording software 'AD')
%
% cdat2waveform(cdat, ...
%  name/value args, any order
%
% one or both of:
%                'upthresh', [], 
%                'downthresh', [],
%
%                'spikewin', 32,
%                'prespikepts', 6,
%                'postignore', 26,
%
% Defaults are designed to mimic Matt Wilson's AD acquisition software.
%
% cdat is a cont struct (as returned by Tom's imcont.m)
%
% -upthresh and downthresh are thresholds in same units as cdat (usually
% mV). If both are specified, we trigger on upward and downward crossings of
% the two values. 'postignore' is enforced across the union of the
% crossings. (E.g. if you trigger on a downward crossing, you have to wait
% 'postignore' samples before you can trigger on either an upward or a
% downward crossing.)
%
% -spikewin is number of waveform points to extract per spike (AD = 32)
%
% -prespikepts is number of points before the threshold crossing to include
% (AD = 6)
%
% -postignore is the number of points after a threshold crossing where we
% ignore other threshold crossings (AD = 26).
%
% Tom Davidson (tjd@alum.mit.edu) 5/25/2006, 1/25/2011

% TODO/wishlist/minutiae
% -different thresholds per chan? (array of thresholds)

% struct cont has correct members, type, dims?
if ~(isstruct(c) && ...
     ndims(c.data)==2);
  error('Invalid data');
end

a = struct('upthresh',[], ...
           'downthresh',[],...
           'spikewin', 32, ...
           'prespikepts', 6, ...
           'postignore', 26);

a = parseArgsLite(varargin,a);

% validate input args:

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

% transpose so I don't have to refactor the old cont2waveform code below :)
data = c.data;

% c is at least a.spikewin samples
if size(data,1) < a.spikewin
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
  mcth = max(data,[],2);
  % When diff of logical == 1, we have a xing from false to true.
  % If signal starts above thresh, not counted as xing.
  upxings = find(diff(mcth > a.upthresh) == 1);
end

% repeat for downward xings
downxings = [];
if ~isempty(a.downthresh),
  mcth = min(data,[],2);
  downxings = find(diff(mcth < a.downthresh) == 1);
end

% sort them both into one list
xings = sort([upxings downxings]);

% logical array of xings to ignore--assume none to start.
% (Faster than resizing xings for each exclusion)
xingsignore = false(length(xings),1);

% don't use points in last postspikepts (to avoid reading off end of
% data)
xingsignore(find(xings < size(data,1)-postspikepts, 1, 'last')+1:end) = true;

% don't use points in first prespikepts (to avoid reading off beginning of data)
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



%%% throw out ignored xings:

xings = xings(~xingsignore);
  

% initialize output struct
nspikes = size(xings,1);
wf = struct('timestamp',[],'waveform',[]);
wf.timestamp = zeros(1,nspikes,class(c.tstart));
wf.waveform = zeros(size(data,2),a.spikewin,nspikes, class(data));

% fill output struct
for j = 1:nspikes ;
  
	% lookup timestamp for this buffer (keep sync'd with AD's clock)
	% (+1 is because timestamp array is 1-indexed)
        % vectorize?
        wf.timestamp(j) = (xings(j)-1)/c.samplerate+c.tstart;
        
	% grab spikewin pts around threshold xing
	wf.waveform(:,:,j) = data(xings(j)-a.prespikepts:xings(j)+postspikepts-1,:);
end
