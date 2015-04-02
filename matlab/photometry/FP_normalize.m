function [c_out c_Regress bls YFit] = FP_normalize(c_Mag, varargin)
% FP_NORMALIZE: normalize a photometry signal using isosbestic channel
% 
% [c_dFF c_Regress bls YFit filt] = FP_normalize(c_Mag, varargin)
%
% REQUIRED INPUTS:
%  c_Mag: 2-channel cont struct, with signal in ch1, isosbestic control in ch2
%
%  Param/Value pairs.
%  'control_LP_F': transition band for LPF to apply to the control signal
%                  before fitting (Hz; no default, try: [0.1 0.2]);
%
% OPTIONAL INPUTS:
%     'norm_type': dF/F denominator: control channel fit (default), constant 
%                  value, or no normalization (see code for discussion): 
%                  ({'fit'}, 'const', 'none')
%    'rig_baseline_V': power in each channel before plugging in the animal
%                  (will be subtracted before normalizing, default [0 0])
% 'dFF_zero_prctile': re-zero dFF signal to this percentile of the final
%                  signal. Empty means do not re-zero (default: 1)
% 'original_timewin': whether to return all data, including invalid regions
%                  due to filter edge effects (default: false)
%
% NB: the beginning and end of the signal may be invalid due to filter edge
% effects, and are excluded by default. Lower frequency filters will result
% in longer excluded regions.
%
% Example:
%
%   c_Mag = contdemodulate(c_FP_Raw,...
%                          'nsignals', 2,'bandwidth_F', [10 15])
%
%   c_dFF = FP_normalize(c_Mag,...
%                        'control_LP_F', [0.1 0.2],...
%                        'norm_type', 'fit');
%
%   c_dFF = contwin(c_dFF,[],'samps_good');
%
% Tom Davidson, Stanford University, Feb 2015. <tjd@alum.mit.edu>


a = struct(...
    'control_LP_F', [],...
    'norm_type', 'fit',...
    'rig_baseline_V', [0 0],...
    'dFF_zero_prctile', 1,...
    'original_timewin', false);

a = parseArgsLite(varargin, a);

if isempty(a.control_LP_F)
    error('Must provide transition band for lowpass filter ''control_LP_F''. Try: [0.1 0.2]');
end

if isempty(a.rig_baseline_V)
    a.rig_baseline_V = [0 0];
end

if numel(a.rig_baseline_V) ~=2,
    error('Must provide 2 values for ''rig_baseline_V'', or leave it empty');
end

switch c_Mag.units
    case {'V'}, 
        c_Mag.units = 'V';
    case 'mV',
        c_Mag.data = c_Mag.data./1000;
        c_Mag.units = 'V';
    otherwise
        error('Unrecognized ''units'' in c_Mag');
end
%% Normalize using 'isosbestic' channel
% Debleaching below works well. Does not make assumptions 
% about dynamics (other than speed, set by LPF), so we can follow unusual 
% bleaching curves.

% Signal to detrend is in chan #1, isosbestic control is in chan 2;
ch_sig = 1;
ch_iso = 2;

% Remove baseline from each signal
c_Mag.data = bsxfun(@minus, c_Mag.data, a.rig_baseline_V([ch_sig ch_iso]));
c_Mag = contdatarange(c_Mag);

% downsample, apply 'very-low-pass' filter, then upsample
fopt_vlp = mkfiltopt('filttype', 'lowpass', 'F', a.control_LP_F, 'name', 'VLP','atten_db', 40);
c_Mag_F_ds = contfilt(c_Mag, 'filtopt', fopt_vlp, 'autoresample', true);

% upsample smoothed control signal to match original
c_Mag_F = continterp(c_Mag_F_ds,...
      'timewin', [c_Mag.tstart c_Mag.tend], ...
      'samplerate', c_Mag.samplerate,...
      'extrapval', 0);
  
% Select channels to regress 
% (Signal to detrend is in chan #1, isosbestic control is in chan 2)
c_Regress = contcombine(...
    contchans(c_Mag,'chans', ch_sig),... 
    contchans(c_Mag_F, 'chans', ch_iso));
c_Regress.chanlabels = {'Y_signal', 'XX_isosbestic_control'};

% Filtering results in edge effects--we don't want to include this data in
% our fit. (contfilt/continterp/contcombine keep track of this
% automatically)
c_Regress = contwin(c_Regress, [], 'samps_good');

% regress FP signal against 405 control to learn coeffs:
% (if numerical problems, try subsampling data)
XX = c_Regress.data(:,ch_iso); % smoothed control channel
Y = c_Regress.data(:,ch_sig); % original signal channel to fit
bls = polyfit(XX,Y,1);

% Model all data (including invalid points not used in the regression)
XX_all = c_Mag_F.data(:,ch_iso); % smoothed control
Y_all = c_Mag.data(:,ch_sig); % original signal

Y_fit_all = bls(1) .* XX_all + bls(2);

% Subtract Y_fit to get the residual 'transients' (in detector units,
% i.e. Volts)
Y_dF_all = Y_all - Y_fit_all;

% 2 options for normalization of deltaF:
%
switch a.norm_type

    case 'const'
        % OPTION 1:
        % constant scaling assumes that GCaMP residual doesn't bleach, so
        % transients don't need to be scaled
        
        % scale deltaF by mean of (valid points in) signal.
        Y_dFF_all = Y_dF_all ./ mean(Y);
        suffix = 'dFF';

    case 'fit'
        % OPTION 2:
        % dividing through by Y_fit assumes that GCaMP bleaching scales with
        % overall tissue bleaching.
        Y_dFF_all = Y_dF_all./Y_fit_all; % deltaF/F for each timepoint
        suffix = 'dFF';

    case 'none'
        % OPTION 3:
        % Don't scale; return dF, in detector units        
        Y_dFF_all = Y_dF_all;
        suffix = 'dF';

        
    otherwise
        error('Unrecognized ''norm_type'' argument');
end

% Optionally shift baselin (to nth percentile of valid data) to make 
% comparable to usual deltaF/F measures
if ~isempty(a.dFF_zero_prctile);
    Y_dFF_valid = Y_dFF_all([c_Mag_F.nbad_start+1:end-c_Mag_F.nbad_end]);
    Y_dFF_all = Y_dFF_all - prctile(Y_dFF_valid,a.dFF_zero_prctile);
end

% save out cont struct:
c_out = mkcdat(...
    'name', [c_Mag.name '_' suffix],...
    'chanlabels', {[c_Mag.chanlabels{ch_sig} '_' suffix]},...
    'samplerate', c_Mag.samplerate,...
    'tstart', c_Mag.tstart,...
    'tend', c_Mag.tend,...
    ... % c_Mag_F signal used to make dFF, has largest invalid regions
    'nbad_start', c_Mag_F.nbad_start,...
    'nbad_end', c_Mag_F.nbad_end,...
    'max_tserr', c_Mag.max_tserr,...
    'data', Y_dFF_all);
c_out = contdatarange(c_out);

if ~a.original_timewin,
    c_out = contwin(c_out, [], 'samps_good');
end

contcheck(c_out);
