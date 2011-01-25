function tth = eeg2tt(eegfiles, ttfile, trodetype, varargin)
% EEG2TT - create a .tt file from an .eeg file (mimics AD)
%   $Id: eeg2tt.m 2213 2009-08-03 19:38:21Z tjd $
%
% eeg2tt(eegfiles, ttfile, 'octrode'/'tetrode1'/'tetrode2', ...
%
%        followed by param/value args, any order
%
%    'upthresh_uv', [],
%    'downthresh_uv, [],
%    'spikewin', 32, 
%    'prespikepts', 6, 
%    'postignore', 26,
%    'ignorechan', [],
%    'chunksize', 20e6
%    'nocheck', false
%    'overwrite', false
% sample usage:
%
% eeg2tt({'f01.etd' 'f01.001' 'f01.002'}, 'f01.tt', 'tetrode2', ...
%  'upthresh_uv', 55, 'downthresh_uv', -80, 'spikewin', 64, 'postignore', 15)
%
% EEG2TT will extract short waveforms from several continuously-recorded
% ('.eeg') data files. Can be useful if you need triggering that is not
% supported by AD (longer spike windows, overlapping spike windows, upward
% and downward-going thresholds).
%
% To precisely mimic the normal function of recordings made using AD, use:
% spikewin of 32 samples, prespikepts of 6 samples, postignore of 26
% samples, and an upward-going threshold (upthresh_uv).
%
% *** arguments:
%
% -eegfiles is cell array of eeg file names, *in order*
% -ttfile is name of new tt tile to create
%
% -upthresh_uv / downthresh_uv are thresholds in microvolts for either
% upward- or downward-going threshold crossings. You can specify both
% (see 'help cont2waveform' for details.)
%
% -spikewin is number of samples to extract per channel per spike (AD=32)
%
% -prespikepts is number of samples before the threshold crossing to include
% in the spikewin (AD=6)
%
% -postignore is the number of samples after a threshold crossing during
% which we ignore other threshold crossings (AD = 26).
%
% -ignorechan is a list of channels to be 'zeroed out' before processing
% for threshold--useful for channels with artifacts. Note this is
% 1-indexed for the electrode being worked on (so 1-4 for either
% tetrode1 or tetrode2)
%
% -chunksize is the size in bytes of the 'chunk' of eegfile to process
% (conserves memory; default: 20e6 = 20MB)
%
% -nocheck: If true, proceed even if the gains are different across input
% files.
%
% -overwrite: If true, proceed even if the output file already exists

% known issues:
% -misses spikes that cross file boundaries (1 every 2GB. Deal with it!)

% wishlist:
% -check valid arguments (eeg filetype...)
% -check that metadata (rate, gain) is same across all files
% -different thresholds / chan? (have to implement in cont2waveform)

% print some debug statements if true
debug = false;

% defaults
args = struct(...
    'upthresh_uv', [], ...
    'downthresh_uv', [], ...
    'spikewin', 32, ...
    'prespikepts', 6, ...
    'postignore', 26, ...
    'ignorechan', [], ...
    'chunksize', 20e6, ...
    'nocheck',false,...
    'overwrite', false);

args = parseArgsLite(varargin,args);

% Electrode type?
if ~ischar(trodetype)
	error ('Electrode type must be string');
end

if ischar(trodetype)
  switch lower(trodetype)
   case 'octrode'
    chans = 1:8;
   case 'tetrode1'
    chans = 1:4;
   case 'tetrode2'
    chans = 5:8;
   case 'stereo1'
    chans = 1:2;
   case 'stereo2'
    chans = 3:4;
   case 'stereo3'
    chans = 5:6;
   case 'stereo4'
    chans = 7:8;
   otherwise
    error(['Unsupported electrode type string:' trodetype]);
  end

elseif isnumeric(trodetype)
    chans = trodetype;

else
  error(['bad ''trodetype'' / chans argument, try: octrode/tetrode1.../stereo1... ' ...
         'or list of channel #s']);
end
    
if ~isempty(args.ignorechan),
  if any(args.ignorechan > length(chans)) || ...
        any(args.ignorechan < 1),
    error(['all ''ignorechan'' values must be between 1 and the number of ' ...
           'channels on the requested electrode)']);
  end
end

if args.overwrite
  writemode = 'overwrite';
else
  writemode = 'write';
end

spikewin = args.spikewin;
prespikepts = args.prespikepts;
postignore = args.postignore;
chunksize = args.chunksize;

if ~iscell(eegfiles)
  eegfiles = {eegfiles};
end
nfiles = length(eegfiles);

gains = zeros(nfiles,length(chans));
[rates recordsizes] = deal(zeros(nfiles,1));

% check that all are eeg files
for j = 1:nfiles;
  try
    fh = mwlopen(eegfiles{j});
  catch 
    error([lasterr sprintf('\n') 'Did you run adextract ?']);
  end
  if ~isa(fh, 'mwleegfile'),
    error(['File is not an .eeg file (class: ' class(fh) ')']);
  end

  fhdr = fh.header;
  
  % check that gains/rates/record sizes same across all files
  
  gtmp = getHeaderGains(fh)';
  gains(j,:) = gtmp(chans);
% $$$   for k = 1:length(chans),
% $$$     gains(j,k) = str2double(fhdr(2).(['channel ' num2str(chans(k)-1) ' ampgain']));
% $$$   end
  rates(j) = str2double(fhdr(2).rate);
  recordsizes(j) = str2double(fhdr(2).buffer_size);
end

if any(diff(gains(:))),
  disp(['Gain mismatch:' sprintf('\n')]);
  disp(gains);
  if ~args.nocheck,    
    error('Gain mismatch, override with ''nocheck'' option');
  else
    disp('Proceeding... (''nocheck'' enabled).');
  end
end

if any(diff(rates)),
  disp(rates);
  error('Rate mismatch');
end

if any(diff(recordsizes)),
  disp(recordsizes);
  error('Record (buffer) size mismatch');
end



% get info from 1st file
fh = mwlopen(eegfiles{1});
fhdr = fh.header;
headerrate = str2double(...
    fhdr(2).rate / fhdr(2).nchannels);
gain = str2double(...
    fhdr(2).(['channel ' num2str(chans(1)-1) ' ampgain']));
recordsize = fh.recordsize;

% adthresh = thresh in uv * 1E-6 V/uv * gain ADVolts/V * 2048/10 ADUnits/ADVolt
upthresh_ad = args.upthresh_uv * 1e-6 * gain * 2048/10;
downthresh_ad = args.downthresh_uv * 1e-6 * gain * 2048/10;

% create empty tt file
tth = mwlcreate(ttfile, 'waveform', ...
                'NChannels',length(chans),...
                'NSamples', spikewin,...
                'mode', writemode);
tth = closeHeader(tth); % enter 'append' mode

wbh = [];

for fileno = 1:nfiles;
  
  fh = mwlopen(eegfiles{fileno});
  
  nrecords = fh.nrecords;

  if chunksize < recordsize;
    warning ('Requested chunk size smaller than single buffer record, using single buffer');
    chunksize = recordsize;
  end
  
  chunks = floor(0:(chunksize/recordsize):nrecords-1);
  if chunks(end) ~= nrecords-1,
    chunks = [chunks nrecords-1];
  end
  
  for chunkno = 1:(length(chunks)-1);

    if(debug); disp(['working on chunk ' num2str(chunkno) ' of ' num2str(length(chunks)-1)]); end

    cont = load(fh,'all',chunks(chunkno):chunks(chunkno+1));

    % select appropriate channels
    cont.data = cont.data(chans,:,:);

    % zero out channels requested to be ignored.
    if ~isempty(args.ignorechan),
      cont.data(:,args.ignorechan) = 0;
    end
    
    % let cont2waveform figure out sampling rate, except in
    % 1-buffer case
    if length(cont.timestamp) < 2,
      rateperchan = headerrate;
    else
      rateperchan = [];
    end

    % ignore first buffer on all but first chunk in a file
    tth = appendData(tth, ...
                     cont2waveform(cont, ...
                                   'upthresh', upthresh_ad, ...
                                   'downthresh', downthresh_ad, ...
                                   'spikewin',spikewin, ...
                                   'prespikepts',prespikepts, ...
                                   'postignore',postignore, ...
                                   'rateperchan',rateperchan, ...
                                   'ignorefirstbuff', (chunkno ~= 1)));

    wbh = mkwaitbar (chunkno/(length(chunks)-1),...
                     fh.filename,fileno,nfiles, wbh);
  
  end

end

% close waitbar window
close (wbh);


%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS

function wbh = mkwaitbar (frac,fname,fileno,nfiles,wbh)

if isempty(wbh),
  wbh = waitbar(0);
end

wbh = waitbar(frac, wbh,...
              ['eeg2tt progress for file ' fname ...
               ' ( file # ' num2str(fileno) '/' num2str(nfiles) ...
               ' )']);

set(wbh,'name','eeg2tt');




