function s = sem(X)
% SEM calculate the standard error of the mean of the samples in X (columnwise)
  
  if any(isnan(X(:))),
    error('X contains NaNs--use ''nansem''');
  end
  
  s = std(X,0,1)./sqrt(size(X,1));