function [data segs_samp] = contsegdatawin(c, t, offset, varargin)
% CONTSEGDATAWIN get windows of data around time points, return as a matrix
%  [data segs_samp] = contsegdatawin(c, t, offset)
%
% Inputs:
%  c - cont struct
%  t - time points
%  offset - start/end time offsets from t to make segs
%
% Outputs: (same as contsegdata)
%  data - data selected from c.data, 1 column per channel
%  segs_samp - m x 2 array of indexes into c.data for each seg

  a = struct(...
      'excludenans', []);
  a = parseArgsLite(varargin,a);

  if numel(offset) ~= 2,
    error ('''offset'' must be a 2-element vector of times');
  end
  
  segs = bsxfun(@plus, t(:),offset(:)');
  
  [data segs_samp] = contsegdata(c, segs, a.excludenans);
