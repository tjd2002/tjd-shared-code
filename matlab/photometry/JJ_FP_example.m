% Sample code to import, demodulate, and normalize raw photometry signals
% recorded with Jeff and Ada's 2-color, 2-frequency rig in the de Lecea lab
%
% Tom Davidson (tjd@alum.mit.edu) Jan 2016


%% set up analysis parameters 
% see function documentation for information on what each param means

% Uncomment one of the blocks below to specify sample data to process

% Example1: short recording with finger over tip of fiber, from JJ via
% email, sampled at 10kHz. Contains raw detector signal, outputs from SRS
% lockins, and copies of the LED modulation signals.
exptname = '012716 10kHz finger test';
path_to_data = '/home/tjd/tmp/';
time_file = '012716_time.mat';
data_file = '012716_data.mat';
data_chanlabels = {'Lockin1', 'Lockin2', 'Det1', 'Ref1X', 'Ref2X'};
signal_labels = {'470nm' '405nm'}; % labels for channels '1' and '2'
carrier_F = [211 531]; % approximate carrier frequencies


%% demodulation parameters (see contdemodulate.m for documentation)
cfg.demod_BW_F = [5 10]; % bandwidth (Hz)
cfg.demod_ripp_db = 0.1; % filter design parameter: bandpass ripple
cfg.demod_atten_db = 60; % filter design parameter: stopband rejection

% normalization parameters (see FP_normalize.m for documenation)
cfg.FPnorm_norm_type = 'fit'; % type of normalization to perform
cfg.FPnorm_control_LP_F = [2 3]; % low-pass filter transition band (Hz)
cfg.FPnorm_dFF_zero_prctile = []; % [] = do not re-zero

% baseline rig fluorescence for each channel with animal not plugged in, 
% in Volts (use [] if you didn't measure this).
cfg.rig_baseline_V = [];


%% Load in data
if ~exist('cache', 'var'), cache = mkcache(); end;

fprintf('Loading\n');
clear data time;
load([path_to_data '/' data_file], 'data');
load([path_to_data '/' time_file], 'time');

c_Raw = imcont('data', data, ...
    'timestamp', time,...
    'name', exptname,...
    'chanlabels', data_chanlabels,...
    'dataunits', 'V');
fprintf('-->Done!\n');

disp(c_Raw);

%%
fprintf('Demodulating raw photometry signal ...\n');
% Demodulate raw detector signal.
[c_Mag, FP_Ref_F, FP_PSDs, cache, c_Raw_out] = ...
    contdemodulate(c_Raw, ...
    'nsignals', 2,...
    'signal_labels', signal_labels,...
    'bandwidth_F', cfg.demod_BW_F,...
    'ripp_db', cfg.demod_ripp_db,...
    'atten_db', cfg.demod_atten_db,...
    'recover_carriers_from_signal', false,...
    'cache', cache);
fprintf('-->Done!\n');
%     'carrier_F', carrier_F,...
%     'carrier_F_winsize', 5,... % What window to look in for peak 

%% ...Plot it
figh_FP = figure(31);
set(figh_FP, 'Position', [0 0 1200 800]);

timewin = [];

sp(1) = subplot(8,1,1:2); 
title({exptname 'Raw detector signal'});
quickplot(contchans(c_Raw_out, 'chanlabels', {'Det1'}),...
    'color', [0.3 0.3 1.0; 1 0.3 1],...
    'subsample', false, ...
    'timewin', timewin);

ylabel({'Signal amplitude' '(V at Detector)'})
xlabel([]);
set(gca, 'TickDir', 'in');
set(gca, 'XTickLabel', '');
yl = ylim;
ylim([0 yl(2)*1.1]); % show y=0

sp(2) = subplot(8,1,3); 
title('Reference 1 signals');
quickplot(contchans(c_Raw_out, 'chanlabels', {'Ref1X' 'Ref1Y'}),...
    'color', [0.3 0.3 1.0],...
    'subsample', false, ...
    'timewin', timewin);
ylabel('Volts to LED')
xlabel([]);
set(gca, 'TickDir', 'in');
set(gca, 'XTickLabel', '');
yl = ylim;
ylim([-0.1*yl(2) yl(2)*1.1]); % show y=0

sp(3) = subplot(8,1,4); 
title('Reference 2 signals');
quickplot(contchans(c_Raw_out, 'chanlabels', {'Ref2X' 'Ref2Y'}),...
    'color', [1 0.3 1],...
    'subsample', false, ...
    'timewin', timewin);
ylabel('Volts to LED')
xlabel([]);
set(gca, 'TickDir', 'in');
set(gca, 'XTickLabel', '');
yl = ylim;
ylim([-0.1*yl(2) yl(2)*1.1]); % show y=0

sp(4) = subplot (8,1,5:6); 
quickplot(contchans(c_Raw, 'chanlabels', {'Lockin1', 'Lockin2'}),...
    'color', [0.3 0.3 1.0; 1 0.3 1],...
    'subsample', false);
title('Output of SRS Lockin amplifiers');
ylabel({'Volts'});
xlabel([]);
set(gca, 'TickDir', 'in');
set(gca, 'XTickLabel', '');


sp(5) = subplot (8,1,7:8); 
quickplot(contchans(c_Mag, 'chanlabels', {'470nm_demod', '405nm_demod'}),...
    'color', [0.3 0.3 1.0; 1 0.3 1],...
    'subsample', false);
title('Output of software demodulation');
ylabel({'Volts'});

linkaxes(sp, 'x')



%% Now attempt some normalization

fprintf('Normalizing photometry signal ...\n');
[c_dFF, c_Regress, bls] = ...
    FP_normalize(c_Mag, ...
    'norm_type', cfg.FPnorm_norm_type ,...
    'control_LP_F', cfg.FPnorm_control_LP_F, ...
...%    'original_timewin', cfg.FPnorm_original_timewin,...
    'rig_baseline_V', cfg.rig_baseline_V,...
    'dFF_zero_prctile', cfg.FPnorm_dFF_zero_prctile);

fprintf('-->Done!\n');

%% ... and Plot it

% scale the raw inputs according to the regression (for plotting)
c_Mag_scaled = contcombine(...
    contchans(c_Regress, 'chans', 1),...
    contfn(contchans(c_Regress, 'chans', 2), 'fn', @(x) x.*bls(1)+bls(2)));

figh_fpnorm = figure(32);
set(figh_fpnorm, 'Position', [0 0 1000 800]);
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
