function opt = mkenvopt(varargin)
% MKFILTOPT makes filter options objects for use by mkfilt
%
% Args:
%  'method', {'hilbert', 'hilbert_complex', 'peaks', 'rms'}
%  'rms_window_t', [] time to average for rms method

opt = struct(...
    'method', 'peaks',...
    'rms_window_t', []);

opt = parseArgsLite(varargin,opt);

if ~any(strcmp(opt.method,{'hilbert', 'hilbert_complex', 'peaks', 'rms'})),
  error('bad ''method''');
end

if strcmp(opt.method, 'rms') && isempty(opt.rms_window_t);
  error('''rms_window_t'' must be provided for ''rms'' method');
end
