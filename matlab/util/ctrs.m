function c = ctrs(arr)
% CTRS return the centers of an array of edges (1D only!)
  
  if ~isvector(arr),
    warning('ctrs treats multi-d arrays as if they were vectors');
  end
  
  c = mean([arr(1:end-1); arr(2:end)]);