function fopts = filtoptdefs(name)
% FILTS_LFP generates common filters for eeg (theta, beta, ripple...)
%
% just calls mkfiltopt

% All frequency values are in Hz.

% make a bunch of bandpass filters

k = 1;
fname{k} = 'ripple';
F{k} = [100 150 250 300];

k = k+1;
fname{k} = 'gamma';
F{k} = [20 30 80 100];

k = k+1;
fname{k} = 'beta';
F{k} = [7 14 28 36];

k = k+1;
fname{k} = 'theta';
F{k} = [6 7 10 12];

for m = 1:k,
  fopts.(fname{m}) = mkfiltopt('name', fname{m},...
                       'filttype', 'bandpass',...
                       'F', F{m}); 
end

m = m+1;
fopts.spw = mkfiltopt('name', 'SPW',...
                      'filttype', 'lowpass',...
                      'F', [30 35]);

for smooth_sd_ms = [1 2 2.5 3 4 5 6 7 8 9 10 12.5 15 17.5 20 22.5 25 ...
                    27.5 30 35 40 45 50 60 70 80 90 100 125 ...
                    150 175 200 225 250 300 350 400 450 500 600 700 800 ...
                    900 1000 1100 1200 1250 1300 1400 1500 1600 1700 1750 ...
                    1800 1900 2000 2250 2.5e3 5e3 10e3]
  m = m+1;
  fname = ['smooth_sd_' num2str(smooth_sd_ms) 'ms'];

  % field names can't have decimal points, replace them with '_'
  fname_valid = fname;
  fname_valid(strfind(fname, '.')) = '_';

  fopts.(fname_valid) = ...
      mkfiltopt('name', ['smooth_sd_' num2str(smooth_sd_ms) 'ms'],...
                'filttype', 'gausswin',...
                'sd_t', smooth_sd_ms ./ 1000);
end

% just return one filter opts, if requested
if exist('name', 'var') && ~isempty(name)
  fopts = fopts.(name);
end