%% Load in data and set up some variables 

% choose a directory containing .ncs files
datadir = 'C:/recordings/OFT/2010-09-16_20-32-04_PD_rat4/';

% load in some LFP data (change filenames to match yours)
cdat_lfp = NlxCSC2cont('CSCDir', datadir, 'CSCFilenames', {'LFP5.ncs' 'LFP6.ncs' 'LFP9.ncs'});

% load in events (e.g. laser pulses)
ev = NlxLoadEV( [datadir 'Events.nev'] );

% select one channel and a 100 second chunk
cdat_lfp5 = contchans(cdat_lfp, 'chanlabel', 'LFP5');
cdat_lfp5_100s = contwin(cdat_lfp5, [24000 24100]);

% same thing as a 1-liner
cdat_lfp5_100s = contwin(contchans(cdat_lfp, 'chanlabel', 'LFP5'), [24000 24100]);


%% Let's filter in beta (15-30Hz)

% create a 'filter options' (fopt) struct
fopt_beta = mkfiltopt('name', 'beta', 'filttype', 'bandpass', 'F', [11 15 30 35]);

% design and apply the specified filter
[cdat_lfp5_100s_Fbeta Fbeta] = contfilt(cdat_lfp5_100s, 'filtopt', fopt_beta);


%% plot filtered and unfiltered vs time

figure;
hold on;
plot(contts(cdat_lfp5_100s),cdat_lfp5_100s.data, 'k');
plot(contts(cdat_lfp5_100s_Fbeta),cdat_lfp5_100s_Fbeta.data, 'r');

%% create binned time data series
%c = imconthist([name/value pairs])