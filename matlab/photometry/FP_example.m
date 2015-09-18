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

tanksdir = ['/Volumes/DATA_EXT2/FP-Protocols-Raw/DataTanks/'];
tankname = 'FP-Vid-RX8-Jul2015_DT1_072115';
blockname = 'Blk-11';
Raw1_chanlabels = {'Det1', 'Ref1X', 'Ref1Y', 'Ref2X', 'Ref2Y', 'ExcitPower'};
signal_labels = {'480nm' '405nm'};

% tanksdir = ['/Volumes/DATA_EXT2/'];
% tankname = 'Kevin';
% blockname = 'FP_workshop_GCaMP_2';
% Raw1_chanlabels = {'Det1', 'Ref1X', 'Ref1Y', 'Ref2X', 'Ref2Y'};
% signal_labels = {'480nm' '405nm'};

% labels for each channel in store 'Raw1'
%Raw1_chanlabels = {'Det1', 'Ref1X', 'Ref1Y', 'Ref2X', 'Ref2Y'}; 

% parameters to contdemodulate.m
cfg.demod_BW_F = [10 15];
cfg.demod_ripp_db = 0.1;
cfg.demod_atten_db = 50;

% parameters to FP_normalize.m
cfg.FPnorm_norm_type = 'fit';
cfg.FPnorm_control_LP_F = [2 3];
cfg.FPnorm_original_timewin = false;
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
    'original_timewin', cfg.FPnorm_original_timewin,...
    'rig_baseline_V', cfg.rig_baseline_V,...
    'dFF_zero_prctile', cfg.FPnorm_dFF_zero_prctile);
fprintf('-->Done!\n');

%% Plot it
subplot 211; quickplot(c_Mag);
subplot 212; quickplot(c_dFF);
