function [lindex time_ind_val] = contlocalmax(cont)
% CONTLOCALMAX find local maxima of each channel of a contdata struct,
%
%   [lindex time_ind_val] = contlocalmax(c)
% 
% Inputs:
%  c - cont struct
%
% Outputs:
%  lindex - logical array index into c.data, true at each peak (see: localmax.m)
%  time_ind_val - cell array of m x 3 matrices, 1 per channel, with one
%      row per peak found, containing the peak time, the index of the peak,
%      and the value at the peak.

% Tom Davidson <tjd@stanford.edu> 2003-2010
  
  [lindex ind_val] = localmax(cont.data);
  
  % prefix the time of each event to the list, so each event is
  % represented on a row as: time, index, value
  if nargout > 1
    for k = 1:length(ind_val),
      if ~isempty(ind_val{k})
        time_ind_val{k} = zeros(size(ind_val{k},1),size(ind_val{k},2)+1);
        time_ind_val{k}(:,1) = ((ind_val{k}(:,1)-1)./cont.samplerate) + cont.tstart;
        time_ind_val{k}(:,[2 3]) = ind_val{k};
      else
        time_ind_val{k} = [];
      end
    end
  end