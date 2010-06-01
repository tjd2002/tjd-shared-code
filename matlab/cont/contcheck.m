function contcheck(c)
% CONTCHECK - private helper function to verify integrity of cont struct

% Tom Davidson <tjd@stanford.edu> 2003-2010

%% Keep code here light--we'd like to call this function a lot

% Check size of c.data
datasz = size(c.data);

if numel(datasz) > 2,
  error([mfilename ':BadDataDims'], 'c.data must be [m x n] array.');
end

nsamps = datasz(1);
nchans = datasz(2);

% Do we have correct # of samples in data array?
nsamps_expected = round(((c.tend-c.tstart)*c.samplerate) + 1);

if nsamps ~= nsamps_expected
  error([mfilename ':BadDataSize'], ...
        ['Wrong number of samples in c.data. ' ...
         'Based on tstart/tend/samplerate, expected %d; found %d.'],...
        nsamps_expected, nsamps);
end


% Do we have correct # of chanvals/chanlabels (either 0, or nchans)?
switch numel(c.chanlabels)
 case 0
  %ok
 case nchans
  %ok
 otherwise
  error([mfilename ':BadNChanlabels'], 'Wrong number of chanlabels in');
end

switch numel(c.chanvals)
 case 0
  %ok
 case nchans
  %ok
 otherwise
  error([mfilename ':BadNChanvals'], 'Wrong number of chanvals.');
end