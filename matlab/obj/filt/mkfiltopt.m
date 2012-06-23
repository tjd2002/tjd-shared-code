function opt = mkfiltopt(varargin)
% MKFILTOPT makes filter options objects for use by mkfilt
%
% Args:
%  'name': string
%  'Fs': sampling frequency (optional, can also be provided when making filter)
%  'filttype':
%  'highpass'/'lowpass'/'bandpass'/'bandstop'/'gausswin'/'rectwin'/'hatwin'/'coeffs'
%
%  'sd_t': length in time of 1 std. dev of the gaussian window (gausswin/hatwin only)
%  'length_t': length in time of gaussian/rect./hat window 
%
%  Below options only for high-/low-/bandpass/bandstop filters
%  'F': design frequencies, in Hz (as for firpmord, firpm)
%       2 freqs for low/highpass (start/end of roll-off)
%       4 freqs for bandpass (starts and ends of roll-off
%       Freqs must be monotonically increasing
%  'atten_db', -50
%  'ripp_db', 1
%  'datatype', 'single' 
%
% All frequency values are in Hz.

opt = struct(...
    'name', [],...
    'Fs', [],...
    'filttype', '',...
    'F', [],...
    'sd_t', [],...
    'length_t',[],...
    'atten_db', -50,...
    'ripp_db', 1,...
    'coeffs', [],...
    'datatype', 'single' ...
    );

opt = parseArgsLite(varargin,opt);

% 'rowize' opt.F
opt.F = opt.F(:)';
opt.filttype = lower(opt.filttype);
opt.atten_db = -abs(opt.atten_db);
opt.ripp_db = abs(opt.ripp_db);

if ~any(strcmp(opt.datatype,{'single', 'double'})),
  error('bad ''datatype''');
end

if any(diff(opt.F) <= 0),
  error('Design frequencies must be increasing');
end

Flen = length(opt.F);
switch opt.filttype,
 case 'coeffs',
  if isempty(opt.coeffs),
    error('''coeffs'' must be provided for filter of type ''coeffs''');
  end
 case 'gausswin'
  if isempty(opt.length_t) && isempty(opt.sd_t);
    error('''length_t'' or ''sd_t'' must be provided for filters of type ''gausswin''');
  end

 case 'rectwin'
  if isempty(opt.length_t);
    error(['''length_t'' must be provided for filters of type ''rectwin'' ' ...
           'or ''hatwin''']);
  end

 case 'hatwin'
  if isempty(opt.sd_t);
    error('''sd_t'' must be provided for filters of type ''hatwin''');
  end

  
 case {'highpass' 'lowpass'}
  if Flen ~= 2,
    error(['highpass and lowpass filters must have 2 design frequencies ' ...
           '(pass/stop or stop/pass)']);
  end
  
 case {'bandpass' 'bandstop'}
  if Flen ~= 4,
    error(['bandpass/bandstop filters must have 4 design frequencies ' ...
           '(stop/pass/pass/stop, or pass/stop/stop/pass)']);
  end
  
 otherwise
  error(['unrecognized filter ''filttype'': try ''highpass'' ''lowpass'' ' ...
         '''bandpass'' ''bandstop'' ''gausswin'' ''rectwin'' ''hatwin''']);
end