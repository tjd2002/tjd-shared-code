function [N bin marg] = histn(X,edges,keepbins,flag)
% HISTN N-dimensional histogram of multivariate data, plus marginal distributions
%
% [N bin marg] = histn (X, edges, [keepbins, 'noedgebin'] )
% 
% -X is the multivariate data, one datapoint per row.
% -edges is a cell array of vectors to use as edges to use for each
% dimension of X. The empty array collects all non-NaN values into one
% bin 
% -keepbins is a cell array of indexes (or logical arrays) into the
% resulting histogram, with bins to keep ([] = keep all)
% -'noedgebin' gets rid of annoying histc behavior of including the last
% provided edge as a bin.
%
% -N: an N-d array of counts
% -bin: the index of the resulting bin for each data point (like histc)
% -marg: cell array of marginal distributions
%
%  (equivalent to histc when called on 1-D data, but doesn't go down columns
%  when only multi-dimensional X and only one edges array provided, like
%  hist/histc does)
  
% test for ascending sort of edges
  
% iterate over each column of X, assigning each value to an index (bin) in the
% output array
  
% accumarray bin into N

  if ~exist('keepbins', 'var')
    keepbins = [];
  end
  
  if ~exist('flag', 'var'),
    noedgebin = 0;
  else
    if isempty(flag),
      noedgebin = 0;
    elseif strcmpi(flag,'noedgebin'),
      noedgebin = 1;
    else
      error('unrecognized flag argument');
    end
  end

  % make compatible with 1-d histc
  if ~iscell(edges),
    edges = {edges};
  end

  nd = length(edges);
  
  if isempty(X),
    bin = zeros(0,nd);
    for dim = 1:nd
      nbins = length(edges{dim}) - (1 * noedgebin);
      marg{dim} = zeros(nbins,1); %#ok
      Nsz(dim) = nbins;
    end
    
    % otherwise zeros makes a square array
    if nd == 1; Nsz(2) = 1; end
    N = zeros(Nsz);
    return
  
  end
  
  if size(X,2) ~= length(edges),
    error(['edges cell array must have same length as number of columns ' ...
           'of X (for 1D data, use X(:))']);
  end
  
  for dim = 1:nd
    if isempty(edges{dim}),
      edges{dim} = [-Inf Inf];
    elseif ~issorted(edges{dim}),
      error('histn edges must be monotonically non-decreasing');
    end
  end
  
  if any(isnan(X(:))),
    warning('histn:NaNwarn', ['Some values of X are NaN--these points will not be included ' ...
             'in output histogram']);
  end
  
  bin = zeros(size(X));
  Nsz = zeros(1,nd);

  % calculate the corresponding indices for each dimension
  for dim = 1:nd
    nedges = length(edges{dim});

    % get size of eventual histogram
    Nsz(dim) = nedges - (1 * noedgebin);

    % convert values in this dim to histbins
    if ~isempty(X),
      [marg{dim} bin(:,dim)] = histc(X(:,dim),edges{dim}); %#ok
    else
      marg{dim} = zeros(nedges,1);
    end
       
    % override annoying histc behavior of including the last edge as a bin
    if noedgebin
      % remove last bin from marginal dist
      marg{dim} = marg{dim}(1:end-1); 
      % convert edge values to 0 (out of range)
      bin(bin(:,dim) == nedges, dim) = 0;
    end
    
  end

  % in 1-D case, N = marginal dist
  if nd ==1,
    N = marg{1};
  
  else % >1-D
    
    % accumulate all values with non-zero indices into the appropriate bin
    %
    % note that accumarray is buggy in R14SP1!, but this bug is not
    % triggered in the case below where we have a 'value' of 1 and use the default
    % accumulation function @sum.
    % http://newsreader.mathworks.com/WebX/.eef34e6?50@632.oZTKbRi72av@

    N = accumarray(bin(all(bin,2),:),1,Nsz);
  end    
  
  if ~isempty(keepbins),
    for dim = 1:nd
      keep{dim} = keepbins{dim}; %#ok
      if isempty(keep{dim}),
        keep{dim} = ':';
      end
    end
    
    N = N(keep{:});
  end