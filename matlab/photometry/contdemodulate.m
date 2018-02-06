function [c_Mag, Ref_F, PSDs, cache, c_Raw] = contdemodulate(c_Raw, varargin)
% CONTDEMODULATE: do amplitude demodulation
%
% [c_Mag, Ref_F, PSDs, cache],...
%      = contdemodulate(c_Raw, <param/value pairs>);
%
% INPUTS:
% c_Raw: cont struct with the following named channels:
%           'Det1': Raw detector signal to be demodulated 
%          'Ref1X': carrier sinusoid for signal 1
%          'Ref1Y': 90-degrees out of phase carrier sinusoid for signal 1
%          'Ref2X': carrier sinusoid for signal 2
%          'Ref2Y': 90-degrees out of phase carrier sinusoid for signal 2
%            ...etc
%
% REQUIRED PARAMETERS:
%       'nsignals': How many signals are we demodulating? (no default)
%    'bandwidth_F': Low-pass transition band, in Hz, for frequencies of 
%                   interest in original signal. Mx2, one row per signal. 
%                   If only one row, use same filter for each signal.
%                   (no default, try [10 15])
%
% OPTIONAL PARAMETERS:
%'detector_chanlabel': Name of detector channel to decode (default: 'Det1')
%        'atten_db': Requested stopband attenuation of filters, in dB 
%                    (default = -50) (NB multiple filters in series)
%         'ripp_db': Requested passband ripple of filters, in dB 
%                    (default = 0.1) (NB multiple filters in series)
%      'LPF_repeat': How many times to apply the low-pass filter after the
%                    product detectors (default = 2)
%'recover_carriers_from_signal': Recover carrier frequencies from input
%                    signal, ignore provided references. (default: false)
%'downsample_outputs': Automatically downsample outputs. (default: true)
%           'cache': object cache containing previously-designed filters,
%                    etc.
%
% OUTPUTS:
%         c_Mag: cont struct, with one channel per demodulated signal
%         Ref_F: array of recovered reference frequencies for each signal
%
% ADDITIONAL DEBUGGING OUTPUTS:
%          PSDs: struct containg spectra of various intermediate signals
%         cache: contains named filters used at various stages
%
% EXAMPLES:
%  [c_Mag RefF] = contdemodulate(c_FP_Raw, ...
%                                'nsignals', 2, ...
%                                'bandwidth_F', [10 15]);
%
% By: Tom Davidson, Stanford University, 2015 (tjd@alum.mit.edu)

% TODO:
% -Handle missing/invalid data (e.g. optostim periods)
%   -NaN/time list? autodetect rails? 
%   1) For each demodulator, censor an integer # of periods of the carrier
%   frequency, then re-insert placeholder values after demodulation.
%     -Possible edge effects due to noise from other channels? 
%     -'Momentum' in LPF before/after censoring
%     -Equivalent to 'pausing' the stream+filtering in TDT code
%   2) Find a LPF that can handle missing values:
%     -KZ, (Kolmogorov-Zurbenko), a series of moving averages. See Wikipedia


    
% -pass out estimated resultant Freqz for all processing.
% -Implement with iFFT instead of bandpass/product detectorr/LPF
% -Generate references from given frequencies (avoid recording RefX/RefY),
%  needed for 'recover_carriers_from_signal', too. (Subtle--can have phase
%  shifts if we get it slightly wrong, but these don't matter for
%  magnitude?)
% -Handle DC as special case/option? Freq 0? Same LPF (code commented out below)
%
% LATER
% -zero-padding/continterp
% -Demodulate multiple detectors (prefilter/choose sampling rates based 
%  on all detectors to avoid aliasing). Maybe implement as a wrapper 
%  function to contdemodulate?
%
% DONE
% -Optionally recover frequencies from raw detector signal?
% -Pass out power spectra of detector, demodulated signals
% -LPF loop
% -Pass out power spectra of detector, demodulated signals


% data integrity check
contcheck(c_Raw);

a = struct (...
    'detector_chanlabel', 'Det1',...
    'nsignals', [],...
    'signal_labels', [],...
    'bandwidth_F', [],...
    'recover_carriers_from_signal', false,...
    'downsample_outputs', true,...
    'atten_db', -50,...
    'ripp_db', 0.1,...
    'LPF_repeat', 2,...
    'cache', mkcache());

a = parseArgsLite(varargin, a);

%% Recover modulation frequencies

% Calculate the input (detector) PSD even if we don't need it below.
c_temp = contchans(c_Raw, 'chanlabels', a.detector_chanlabel);
[PSDs.Pxx_Det, PSDs.F_Det] = periodogram(c_temp.data,blackman(size(c_temp.data,1)),[],c_temp.samplerate);
[PSDs.Pxx_Det_d, PSDs.F_Det_d] = periodogram(detrend(c_temp.data, 'constant'),blackman(size(c_temp.data,1)),[],c_temp.samplerate);

if a.recover_carriers_from_signal, 
    
    % Recover nsignals frequencies with most power from the detector signal
    [~,Ppks_val] = localmax(PSDs.Pxx_Det);
    Ppks_rank = sortrows(Ppks_val{1},-2);
    Ref_F = PSDs.F_Det(Ppks_rank(1:a.nsignals,1));
    
    % Discard unneeded channels from c_Raw:
    c_Raw = contchans(c_Raw, 'chanlabels', a.detector_chanlabel);

    % Have yet to write code to synthesize reference frequencies below
    disp('Demodulation without reference signals not supported yet, returning recovered ref frequencies only.');
    c_Mag = []; cache = [];
    return

else % Recover frequencies from the individual modulation channels
    
    % Ensure we have all the references we'll need
    synthY = false(a.nsignals,1);
    for j = 1:a.nsignals;
        RefXstr{j} = ['Ref' num2str(j) 'X'];
        RefYstr{j} = ['Ref' num2str(j) 'Y'];
        if sum(strcmp(RefXstr{j}, c_Raw.chanlabels)) ~=1,
            error('Must provide reference channels named: ''%s''', RefXstr{j});
        end
        if sum(strcmp(RefYstr{j}, c_Raw.chanlabels)) ~=1,
            synthY(j) = true;
            warning('Quadrature ref signal ''%s'' not provided, will be synthesized', RefYstr{j});
            % error('Must provide reference channel named: ''%s''', RefYstr{j});
        end
    end
    
    % Find peak of FFT in each RefX signal    
    for j = 1:a.nsignals;
        chanidx = chansfromlabels(c_Raw, RefXstr{j});
        [PSDs.Pxx_Ref(:,j), PSDs.F_Ref(:,j)] = ...
            periodogram(detrend(c_Raw.data(:,chanidx)),blackman(size(c_temp.data,1)),[],c_Raw.samplerate);
        [~,maxidx] = max(PSDs.Pxx_Ref(:,j));
        
        Ref_F(j) = PSDs.F_Ref(maxidx); %#ok
    end
    
    % Synthesize missing quadrature (Y) channels if needed
    k = 0;
    if any(synthY),
        for j = 1:a.nsignals,
            if synthY(j),
                
                % we want to delay the reference by 90 deg (1/4 cycle)
                delay = 1/4 * 1/Ref_F(j);
                
                % create RefY as a delayed copy of RefX (by altering tstart/tend)
                k = k+1;
                c_RefjY(k) = contchans(c_Raw, 'chanlabels', RefXstr{j});
                c_RefjY(k).chanlabels = RefYstr(j);
                c_RefjY(k).tstart = c_RefjY(k).tstart+delay;
                c_RefjY(k).tend = c_RefjY(k).tend+delay;
            end
        end
        
        % combine our new reference signals with the originals (contcombine
        % deals with interpolating the delayed signals, and cropping to the
        % valid, overlapping region)
        c_Raw = contcombine(c_Raw, c_RefjY, 'match_first', true);
    end
    
    % Discard unneeded channels from c_Raw:
    c_Raw = contchans(c_Raw, 'chanlabels', {a.detector_chanlabel RefXstr{:} RefYstr{:}});
    
end


%% Analysis parameters setup

% Use same LPF for all signals if only one provided
if size(a.bandwidth_F,1)==1,
    a.bandwidth_F = repmat(a.bandwidth_F,a.nsignals,1);
end
if size(a.bandwidth_F,1)~=a.nsignals,
    error('Must provide either one ''bandwidth_F'', or one per demodulation frequency');
end

for j = 1:a.nsignals;
    % Low-pass filter design parameters (transition band, in Hz)
    fopt_LPF(j) = mkfiltopt(...
        'name', sprintf('LPF%d', j),...
        'filttype', 'lowpass', ...
        'F', a.bandwidth_F(j,:),...
        'atten_db', a.atten_db,...
        'ripp_db', a.ripp_db);

    % Design bandpass filters for each signal to select modulated
    % signal+sidebands (should these be wider?)
    fopt_prefilter(j) = mkfiltopt(... 
        'name', sprintf('BPF%d', j),...
        'filttype', 'bandpass', ...
        'F', Ref_F(j) + [-fliplr(a.bandwidth_F(j,:)) a.bandwidth_F(j,:)],...
        'atten_db', a.atten_db,...
        'ripp_db', a.ripp_db); %#ok
end

% Downsample raw detector channel as needed (to ~6x the max ref frequency +
% upper sideband, so still oversampled). (Was 3x, but forgot about 2f
% rectification)

res_f = double(6*(max(Ref_F)+max(a.bandwidth_F(:,2))) ./ c_Raw.samplerate);

if res_f<c_Raw.samplerate;
c_Raw_ds = contresamp(...
    c_Raw, ...
    'resample', res_f);
else
    c_Raw_ds = c_Raw;
end

% % Skip downsampling for now, until we understand performance of demodulator
% % better
% c_Raw_ds = c_Raw;

c_Det = contchans(c_Raw_ds, 'chanlabels', a.detector_chanlabel);

% %% If this is an expt with no modulation, just filter the raw detector signal

% This should be behind a switch for expts with no modulation?
% Filter the unmodulated detector signal, just for comparison
for j = 1:a.nsignals,
    c_Det_f(j) = c_Det; % initialize
    for k = 1:a.LPF_repeat,
        c_Det_f(j) = contfilt(c_Det_f(j), 'filtopt', fopt_LPF(j), 'autoresample', true);
    end
end

% Hack to save this for output
PSDs.c_Det_f = c_Det_f;

%% Recover magnitude by demodulating.

% Loop over all requested demodulation frequencies
for j = 1:a.nsignals,
    
    % Prefilter with a bandpass before the product detectors to get
    % rid of large signal from other modulated signal.
    % (Don't resample so we can multiply easily with the references)
    [c_Det_pf(j), prefilt(j)] = contfilt(c_Det, ...
        'filtopt', fopt_prefilter(j), ...
        'autoresample', false,...
        'cache', a.cache); %#ok
    a.cache = mkcache('cache', a.cache, 'add_obj', prefilt(j));
    
    % In-phase product detector. pseudocode: X = Det.*RefX (-> LPF)xN;
    % (see SUBFUNCTIONS, below);
    [c_X(j), a.cache] = subf_demodulate(...
        c_Det_pf(j),... % prefiltered detector signal
        contchans(c_Raw_ds,'chanlabels',['Ref' num2str(j) 'X']),...  % RefX
        fopt_LPF(j),... % Lowpass filter specs
        a.LPF_repeat,...
        a.cache); % object cache

    % In-phase product detector. pseudocode: Y = Det.*RefY (-> LPF)xN;
    [c_Y(j), a.cache] = subf_demodulate(...
        c_Det_pf(j),... % prefiltered detector signal
        contchans(c_Raw_ds,'chanlabels',['Ref' num2str(j) 'Y']),... % RefY
        fopt_LPF(j),... % Lowpass filter specs
        a.LPF_repeat,...
        a.cache); % object cache
    
    % Add the estimates of X & Y in quadrature to recover magnitude R
    c_Mag(j) = contoperator(c_X(j), c_Y(j),...
        'op', @hypot); %#ok
    
    fprintf('Done for channel %d, Ref Freq %0.2f Hz\n\n', j, Ref_F(j));
end

c_Mag = contcombine(c_Mag(1), c_Mag(2:end), 'name', [c_Raw.name '_Demod']);

if ~a.downsample_outputs
    c_Mag = continterp(c_Mag, 'samplerate', c_Raw.samplerate);
end

% debug: upsample to original input frequency and take PSD:
c_temp = continterp(c_Mag, 'samplerate', c_Raw.samplerate);
for j = 1:a.nsignals
    [PSDs.Pxx_Mag(:,j), PSDs.F_Mag(:,j)] = periodogram(c_temp.data(:,j),blackman(size(c_temp.data,1)),[],c_temp.samplerate);
end

if ~isempty(a.signal_labels);
    for j = 1:a.nsignals
        c_Mag.chanlabels{j} = [a.signal_labels{j} '_demod'];
    end
else
    for j = 1:a.nsignals
        c_Mag.chanlabels{j} = ['demod_' a.detector_chanlabel '_' sprintf('%0.1f', Ref_F(j)) 'Hz'];
    end
end

cache = a.cache; % for output

% data integrity check
contcheck(c_Mag);

%% SUBFUNCTIONS
function [c_Demod, cache] = subf_demodulate(c_Det_pf, c_Ref, fopt_LPF, LPF_repeat, cache)

% Take product of signal and ref
c_Demod = contoperator(c_Det_pf, c_Ref, ...
    'op', @times, 'infix', '*');

% Apply the post-demodulation LPF as many times as requested
for k = 1:LPF_repeat
    [c_Demod, lpfilt] = contfilt(c_Demod, ...
        'filtopt', fopt_LPF, ...
        'autoresample', true, ...
        'cache', cache);
    cache = mkcache('cache', cache, 'add_obj', lpfilt);
end

  
