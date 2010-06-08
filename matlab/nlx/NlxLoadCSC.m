function [cs NumberValid] = NlxLoadCSC(varargin)
% NlxLoadCSC: gets raw data and headers from .ncs files
%
% Wrapper for Nlx2MatCSC. Loads data, parses headers, converts units 
%
% Tom Davidson <tjd@stanford.edu> 4/2010

% NOTES:
%  -Samples are stored as 16-bit signed integers in the .ncs file, (per
%   www.neuralynx.com/static/software/NeuralynxDataFileFormats.pdf ) but
%   Nlx2Mat returns double-precision floats. Single-precision floats (i.e.
%   the 'single' datatype in Matlab) are sufficient, and use half the
%   memory. Default is the convert to 'single'; this behavior can be
%   controlled with the 'Datatype' argument.
%
%  -The Digital Cheetah acquisition system digitally filters the incoming
%   LFP signals. This introduces a fixed delay of about 1/2 the filter
%   length that is not accounted for when writing samples and timestamps to
%   disk. A single event recorded on one of the digital input channels
%   (which are not filtered), and an LFP or spike channel with different
%   filter settings, will appear to occur at different times. If the
%   'CorrectFilterDelay' option is set to true (the default), we adjust the
%   output timestamps to account for the delay. Cheetah filters before
%   downsampling, so the filtering always takes place at the native
%   sampling rate of 32556Hz. Filters used are (reportedly) linear phase,
%   so the delay introduced is (ntaps - 1) / (2 * 32556). The number of
%   taps in each (low and high) filter is recorded in the header of the
%   .ncs or .ntt files. 
% (See http://www.neuralynx.com/static/newsletters/6.htm for more details).
%
%  -Timestamps are stored as unsigned 64-bit integers (uint64) in the
%   original file, but there's no reason to convert from 'double', since the
%   memory requirements are the same, and the floating point precision
%   (given by 'eps) is << 1 even for very large timestamp values
%   (corresponding to years-long recordings).
%
%  -DON'T TRUST the sampling frequency reported by by Nlx2Mat. Instead, we
%   calculate it from the timestamps given at the beginning of each block of
%   samples (cs.infofS_actual). The value returned by Nlx2Mat for each
%   512-sample record seems to be the sampling rate that was requested
%   during recording, and appears to be an integer. This is useful for
%   detecting a change in settings during recording (as below), but it is
%   not the true sampling frequency. The '-SamplingFrequency' value written
%   in the file header is given to 10 significant digits, but is also wrong
%   (Per MRW, it is derived from the 32556Hz base sampling rate of the
%   cheetah system, and the requested 'interleave'). E.g. in one file, the
%   per-block fS was 1550, the fS in the header was 1550.285714
%   (i.e. 32556/21), and the true fS (calculated from timestamps across the
%   whole file was 1550.099206. So it's best not to use the reported values for
%   synthesizing intermediate timestamps (and definitely not for
%   extrapolating times based on number of elapsed samples), to avoid
%   accumulating small timing errors. Note that imcont.m handles this
%   problem for you; ask Tom for details.
%
%  -If recording or acquisition is stopped and restarted, it becomes harder
%   to calculate the exact sampling frequency. Our approach is to assume that
%   the value reported in the header is close to the true value, and to
%   calculate the true sampling frequency from those chunks of the file that
%   appear to be sampled at near (within 1%) of that frequency. This excludes
%   gaps in recording. It might also be possible to parse the Cheetah
%   logfile for start/stop times.
%
%  -For now, if the user requests a 'TimeWin' that spans a gap in
%   recording, we pad with the value given in the 'PadValue' argument
%   (default NaN). 

% TODO:
%  -missing timestamps cleanup
%    -if requested start/end time is in a gap?
%    -if requested timewin is completely within a gap? use fS_header?
%  -handle gaps without padding (return cell arrays of samps/times?)
%  -import in chunks to save memory? (use indexes rather than timestamps?)
%  -ensure no dependencies to tjd code libraries (subfunctions?)
%
% DONE:
%  -support time ranges (use -Inf/Inf for start/end)
%  -handle gaps in recording (pad)

%% check dependencies

CheckForNlx2Mat('CSC');

%% constants

% nominal native sampling frequency of Digital Cheetah
CheetahfS = 32556;

% Cheetah timestamps per second
ts_per_second = 1e6;

% Size of Neuralynx header, in bytes, per
% http://www.neuralynx.com/static/software/NeuralynxDataFileFormats.pdf
hdr_bytes = 16384; 


%% set up function input arguments
p = inputParser;
p.addRequired('Filename', @sub_isfile); 
p.addParamValue('TimeWin', [], @isnumeric); % in TimeUnits!
                                                    % (usually seconds)
p.addParamValue('Datatype', 'single', @(s)sub_isinlist(s,{'single','double','int16', 'int32'}));
p.addParamValue('DataUnits', 'mV', @(s)sub_isinlist(s,{'ADbits','mV','V'}));
p.addParamValue('TimeUnits', 'seconds', @(s)sub_isinlist(s,{'seconds', 'microseconds'}));
p.addParamValue('UnwrapBuffers', true, @islogical);
% if buffers are missing, should we pad when unwrapping?, and with what value?
p.addParamValue('PadMissingSamples', true, @islogical); 
p.addParamValue('PadValue', NaN, @isnumeric); 
p.addParamValue('CorrectFilterDelay', true, @islogical);
p.addParamValue('IgnoreInvalid', false, @islogical);
p.addParamValue('debug', false, @islogical);

% parse arguments to arglist
p.parse(varargin{:});
a = p.Results;

if ~strcmp(a.DataUnits, 'ADbits') && strncmp('int', a.Datatype, 3),
    error ('Integer datatypes are only valid for ''ADbits'' data units');
end

if ~a.UnwrapBuffers && a.PadMissingSamples,
    warning ('Can only pad missing samples when unwrapping buffers');
end

%% handle requested times
% Get the timestamp of the first and last data buffer in the file, and then
% 'crop' the requested TimeWin so that it doesn't go beyond these limits

% get the total number of buffers in the file
hdr = NlxParseHeader(Nlx2MatCSC(a.Filename, [0 0 0 0 0], true, 1));
finfo = dir(a.Filename);
nrecords = (finfo.bytes - hdr_bytes) ./ hdr.RecordSize;

% get just the first/last buffer timestamp (note zero-indexed)
[buff_tsrange] = Nlx2MatCSC(a.Filename, [1 0 0 0 0], false, 3, [0 nrecords-1]);

% Empty timewin means select all times
if isempty(a.TimeWin),
  a.TimeWin = [-Inf Inf];
end

% Get requested time range in timestamp units
switch a.TimeUnits
 case 'seconds'
  tsrange = a.TimeWin*ts_per_second;
 case 'microsecodns',
  tsrange = a.TimeWin/1e6*ts_per_second;
end

% now handle various problems with the requested times

if tsrange(2) < tsrange(1),
  error('TimeWin end time must not be before start time');
end

if tsrange(1) > buff_tsrange(2) || tsrange(2) < buff_tsrange(1)
  error(sprintf('Requested TimeWin [%0.3f %0.3f] out of range [%0.3f %0.3f].', ... 
                tsrange, buff_tsrange)); %#ok
end

% crop to actual timestamp ranges (Nlx2MatCSC v4 for windows barfs on -Inf/Inf)
if tsrange(1) < buff_tsrange(1),
  tsrange(1) = buff_tsrange(1);
end

if tsrange(2) > buff_tsrange(2),
  tsrange(2) = buff_tsrange(2);
end
  

%% set up extraction parameters and run extraction
FieldSelection = [1 1 1 1 1]; % [timestamps chan fS numvalid samples]
ExtractHeader = 1;
ExtractMode = 4; % 1-all, 2-range, 3-list, 4-timestamp range, 5-ts list
ModeArray = tsrange; % blank for mode 1, 2 elements for range (mode 2 or 4);
                     % n elements for list (mode 3 or 5). indexes (modes 2 and
                     % 3) are zero-indexed

[cs.bufftimes cs.info.chan cs.info.fS_fromfile NumberValid cs.samples cs.info.rawheader] = ...
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

% convert from 1 timestamp/buffer to 1 timestamp/sample
if a.UnwrapBuffers
    [bufflen nbuffs] = size(cs.samples);
    
    % do some heuristic stuff to see if there are missing times in the file due to
    % starting/stopping recording, then try to find the 'true' fS anyways.
    
    % expected interval between buffer start times
    binterval_fromfile = (1 ./ cs.info.fS_fromheader) .* bufflen;
    
    % observed intervals
    bintervals = diff(cs.bufftimes);
    
    % expect there to be less than 1% error in fS reported in header
    fstol = 0.01;

    % no buff intervals should be shorter than expected
    if any(bintervals < binterval_fromfile * 1-fstol);
      error(['Data buffers shorter than expected -- could be due to ' ...
             'significant misreporting of sampling rate in header.']);
    end
    
    % Ff buff intervals are longer than expected, we have gaps
    gapi = bintervals > binterval_fromfile * 1+fstol;

    if any(gapi)
      disp(['Found gaps at (s) : ' num2str(cs.bufftimes([gapi false]))]);
      disp(['Gap durations (s) : ' num2str(bintervals([gapi false]))]);
    end

    % calculate actual sampling rate from *typical* buff interval (i.e. mean of
    % all non-gap intervals). 
    binterval_typ = mean(bintervals(~gapi)); % typical bufftime interval
    sinterval_typ = binterval_typ ./ bufflen; % typical sampling interval
    bufflen_t = bufflen*sinterval_typ; % duration of typical buffer

    % this is our best guess at 'actual' sampling rate
    cs.info.fS_actual = 1./sinterval_typ;
    
    %% extrapolate sample times using bufftimes and fs_actual

    % get recording gaps between buffers (0, not sinterval_typ, for non-gap case)
    bgaps = bintervals-bufflen_t;
    
    % preallocate unwrapped result vector
    if a.PadMissingSamples
      padsamps = ceil(sum(bgaps(gapi))./sinterval_typ);
    else
      padsamps = 0;
    end
    newts = NaN(numel(cs.samples)+padsamps,1);
    newsamps = newts;
    
    if a.debug,
      disp(['initial newts sz: ' num2str(size(newts))]);
    end

    % precompute time increments for a normal buffer
    time_inc = (0:bufflen-1)' .* 1/cs.info.fS_actual; 

    % gapi is index into gaps, make it an index into bufftimes
    gapi = [gapi false];
    
    j = 1; % index into new timestamp result vector
    offset = 0; % accumulator for sampling phase offsets (in samples)
    for k = 1:nbuffs
      newts(j:(j+bufflen-1)) = cs.bufftimes(k) + time_inc;
      newsamps(j:(j+bufflen-1)) = cs.samples(:,k);
      j = j+bufflen;
      % pad the gaps, if requested 
      if a.PadMissingSamples && gapi(k) 

        % Phase of sampling may not be preserved across gaps; we keep track of these
        % phase differences to avoid accumulating timestamp errors greater
        % than 1 sample across multiple gaps. This preserves uniform
        % sampling across entire file (useful for, e.g. imcont) at the
        % expense of uniform sampling across gaps. (During gaps in recording
        % Neuralynx seems to keep the sampling clock going, and restarts
        % recording at same phase; it would take a long gap to accumulate
        % >0.5 sample worth of error; not sure about gaps in acquisition)
        npad = round(bgaps(k)./sinterval_typ);
        
        % calculate offset *in samples*, not seconds
        offset = offset + (bgaps(k)-(npad.*sinterval_typ));
        if a.debug,
          disp(['offset: ' num2str(offset)]);
        end
        if offset >= 1,
          npad = npad+1; % add an 'extra' sample to this gap
          offset = offset-1;
          if a.debug,
            disp('extrapad');
          end
        end
        
        % linspace between last timestamp and start of next buffer
        newts(j-1:(j+npad)) = linspace(newts(j-1),cs.bufftimes(k+1),npad+2);
        newsamps(j:(j+npad-1)) = a.PadValue;

        % newts(j+npad) is the t at start of next buffer, will get overwritten:
        j = j+npad; 

% $$$ %We don't do it this way because if we were to add in an 'extra'
% $$$ %padding ts, we'd get non-monotonic timestamps
% $$$         newts(j:(j+npad-1)) = newts(j-1) + (1:npad) .* sinterval_typ;
% $$$         j = j+npad; 
      
      end
    end

    if isnan(newts(end)) % we preallocated 1 sample too many
      newts(end) = [];
      newsamps(end) = [];
    end
    
    if a.debug,
      disp(['final newts sz: ' num2str(size(newts))]);
      disp(['#NaNs in newts: ' num2str(sum(isnan(newts)))]);
      disp(['issorted newts? ' num2str(issorted(newts))]);
    end
    
    % write times to 'samptimes' and unwrap samples into a vector
    if ~a.debug
      cs = rmfield(cs, 'bufftimes');
    end
    cs.samptimes = newts;
    if a.debug
      cs.oldsamps = cs.samples;
    end
    cs.samples = newsamps;
end


% save out arguments to this function
cs.info.loadargs = a;

% order struct fields alphabetically
cs = orderfields(cs);
cs.info = orderfields(cs.info);


%% Subfunctions
function tf = sub_isfile(s)
if exist(s)~=2
    error('Must be a valid file.');
end
tf = true;

function tf = sub_isinlist(s, validlist)
if ~any(strcmpi(s, validlist))
    error(['Value ''' s ''' not in valid list ' cell2str(validlist)]);
end
tf = true;