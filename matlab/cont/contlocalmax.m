function [lindex time_ind_val] = contlocalmax(cont)
% CONTLOCALMAX apply localmax to each channel of a contdata struct,
% return times of peaks
  
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