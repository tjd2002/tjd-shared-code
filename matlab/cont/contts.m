function [ts] = contts(c)
% CONTTS return timestamps corresponding to data in c
  ts = linspace(c.tstart, c.tend, size(c.data,1));
