% Sample code to import, demodulate, and normalize raw photometry signals
% recorded with a 2-color, 2-frequency rig.
%
% Tom Davidson (tjd@alum.mit.edu) April 2015


%% set up analysis parameters 
% see function documentation for information on what each param means

% set up directory names to find your data
% tanksdir = ['/users/tjd/Data/ChatRat-FP-2015-Raw/DataTanks/'];
% tankname = 'NNLFP-FP-Feb2015_DT2_021315';
% blockname = 'Blk-5';
% Raw1_chanlabels = {'Det1', 'Ref1X', 'Ref1Y', 'Ref2X', 'Ref2Y', 'X1', 'Y1'};

% Edit this path so it points to the 'Example Data' folder from the
% Fiber Photometry Workshop shared folder on 
path_to_Example_Data = '/Users/tjd/Google Drive/Fiber Photometry Workshop/Example Data/';

tanksdir = [path_to_Example_Data '/DataTanks/'];
if ~exist(tanksdir, 'dir'),
    error('Sample data folder not found at: %s', [path_to_Example_Data tanksdir])
end

%% Uncomment one of the blocks below to specify sample data to process

% Example1: ~10min of data from a GCaMP6f-expressing mouse
exptname = 'GCaMP6f-expressing mouse, 20 days post-injection, open field';
tankname = 'FP-Vid-RX8-Jul2015_DT1_072115';
blockname = 'Blk-13';
Raw1_chanlabels = {'Det1', 'Ref1X', 'Ref1Y', 'Ref2X', 'Ref2Y', 'ExcitPower'};
signal_labels = {'480nm' '405nm'};

% % Example2: ~5min of data from a GFP control mouse
% exptname = 'GFP-expressing mouse, 20 days post-injection, open field';
% tankname = 'FP-Vid-RX8-Jul2015_DT1_072115';
% blockname = 'Blk-11';
% Raw1_chanlabels = {'Det1', 'Ref1X', 'Ref1Y', 'Ref2X', 'Ref2Y', 'ExcitPower'};
% signal_labels = {'480nm' '405nm'};


%% demodulation parameters (see contdemodulate.m for documentation)
cfg.demod_BW_F = [10 15]; % bandwidth (Hz)
cfg.demod_ripp_db = 0.1; % filter design parameter: bandpass ripple
cfg.demod_atten_db = 50; % filter design parameter: stopband rejection

% normalization parameters (see FP_normalize.m for documenation)
cfg.FPnorm_norm_type = 'fit'; % type of normalization to perform
cfg.FPnorm_control_LP_F = [2 3]; % low-pass filter transition band (Hz)
cfg.FPnorm_dFF_zero_prctile = []; % [] = do not re-zero

% baseline rig fluorescence for each channel with animal not plugged in, 
% in Volts (use [] if you didn't measure this).
cfg.rig_baseline_V = [];


%% Load in raw FP signals (detector output, carrier frequencies)
if ~exist('cache', 'var'), cache = mkcache(); end;

fprintf('Loading store: ''Raw1''\n');
S = TDT_Import(tanksdir, tankname, blockname, 'Raw1');
c_FP_Raw = imcont('tdtwave', S, ...
    'chans', [],...
    'chanlabels', Raw1_chanlabels,...
    'dataunits', 'V');
fprintf('-->Done!\n');

%%
fprintf('Demodulating raw photometry signal ...\n');
% Demodulate raw detector signal.
[c_Mag, FP_Ref_F, FP_PSDs, cache] = ...
    contdemodulate(c_FP_Raw, ...
    'nsignals', 2,...
    'signal_labels', signal_labels,...
    'bandwidth_F', cfg.demod_BW_F,...
    'ripp_db', cfg.demod_ripp_db,...
    'atten_db', cfg.demod_atten_db,...
    'cache', cache);
fprintf('-->Done!\n');

fprintf('Normalizing photometry signal ...\n');
[c_dFF, c_Regress, bls] = ...
    FP_normalize(c_Mag, ...
    'norm_type', cfg.FPnorm_norm_type ,...
    'control_LP_F', cfg.FPnorm_control_LP_F, ...
    'rig_baseline_V', cfg.rig_baseline_V,...
    'dFF_zero_prctile', cfg.FPnorm_dFF_zero_prctile);

fprintf('-->Done!\n');

%% Plot it

% scale the raw inputs according to the regression (for plotting)
c_Mag_scaled = contcombine(...
    contchans(c_Regress, 'chans', 1),...
    contfn(contchans(c_Regress, 'chans', 2), 'fn', @(x) x.*bls(1)+bls(2)));

figure('Position', [0 0 900 500]);
sp(1) = subplot(3,1,1); 
title(exptname);
quickplot(contfn(c_Mag, 'fn', @(x)x*1000),...
    'color', [0.3 0.3 1.0; 1 0.3 1],...
    'subsample', false);
ylabel({'Fluorescence amplitude' '(mV at Detector)'})
yl = ylim;
ylim([0 yl(2)*1.1]); % show y=0
lh(1) = legend;
lh(1).Location = 'SouthWest';
lh(1).Interpreter = 'none';

sp(2) = subplot (3,1,2); quickplot(contfn(c_Mag_scaled, 'fn', @(x)x*1000),...
    'color', [0.3 0.3 1.0; 1 0.3 1]),...
    'subsample', false;
title('Best fit of smoothed isosbestic control to signal channel');
ylabel({'Fluorescence amplitude'  '(scaled)'});
lh(2) = legend;
lh(2).Location = 'NorthEast';
lh(2).Interpreter = 'none';

sp(3) = subplot (3,1,3); quickplot(c_dFF, ...
    'subsample', false, ...
    'color', [0.3 0.6 0.3]); 
ylim([-0.1 0.25]) % common dF/F range
title('Normalized signal');
ylabel('dF/F');
linkaxes(sp, 'x')


%% Plot 'under the hood'

% scale the raw inputs according to the regression (for plotting)
c_Mag_scaled = contcombine(...
    contchans(c_Regress, 'chans', 1),...
    contfn(contchans(c_Regress, 'chans', 2), 'fn', @(x) x.*bls(1)+bls(2)));

figure('Position', [0 0 900 500]);
sp(1) = subplot(3,1,1); 
title(exptname);
quickplot(contfn(c_Mag, 'fn', @(x)x*1000),...
    'color', [0.3 0.3 1.0; 1 0.3 1],...
    'subsample', false);
ylabel({'Fluorescence amplitude' '(mV at Detector)'})
yl = ylim;
ylim([0 yl(2)*1.1]); % show y=0
lh(1) = legend;
lh(1).Location = 'SouthWest';
lh(1).Interpreter = 'none';

sp(2) = subplot (3,1,2); 
quickplot(contfn(contchans(c_FP_Raw, 'chans', [1 2 4]), 'fn', @(x)x*1000), ...
        'subsample', false);
title('Raw references and detector signal');
ylabel({'mV'});
lh(2) = legend;
lh(2).Location = 'NorthEast';
lh(2).Interpreter = 'none';

sp(3) = subplot (3,1,3); quickplot(c_dFF, ...
    'subsample', false, ...
    'color', [0.3 0.6 0.3]); 
ylim([-0.1 0.25]) % common dF/F range
title('Normalized signal');
ylabel('dF/F');
linkaxes(sp, 'x')

