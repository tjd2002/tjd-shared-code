function [lindex varargout] = localmax(x)
% LOCALMAX find indexes of local maxima in an array
%from http://www.mit.edu/~pwb/cssm/matlab-faq.html
%
% **operates column-wise on matrices
%
%  [lindex values] = localmax(x)
%
%  lindex: a logical index into x which is true at local maxima of x
%  values: cell array of index/value pairs at the local maxima for each
%  column (same # of rows as sum(lindex) for that column)
  
%  lindex = diff( sign( diff([0; x(:); 0]) ) ) < 0;
%
%  unfortunately this doesn't deal well with plateaus (either flat tops, or
%  intermediate repeating values). We want to : 1) If it's a true peak,
%  report the beginning of the plateau. 2) If it's not a true peak, don't
%  report it.
%  
%  We implement this by doing the peak calculation with repeated values
%  omitted (but have to be careful about getting indexes right for output)
  
 
  if ndims(x) > 2, 
    error('localmax only supports up to 2-dimensional arrays');
  end
  
  if size(x,1) == 1, 
    warning('localmax operates columnwise');
  end
  
  [nrows ncols] = size(x); %#ok

  % init logical index of peaks
  lindex = false(size(x));

  for k = 1:ncols
    % have to operate in a for loop b/c we resize the array to delete plateaus
    
    xk = x(:,k);
    
    % get indexes of points in x  making plateaus (i.e. second and later points
    % in a plateau)
    xkplateaui = sign(diff([0; xk(:)])) == 0;
    
    if any(xkplateaui) % possibly save memory by avoiding big assignment
                      % delete these points from x
      xk(xkplateaui) = [];
    end
    
    % calculate peaks, insert them into correct position in lindex
    if any(~xkplateaui)
      lindex(~xkplateaui,k) = diff( sign( diff( [0; xk(:); 0]))) < 0;
    end
    
  end

  % not sure how this becomes double, but we want logical array for indexing
  lindex = logical(lindex);
  
  % only compute this if requested
  if nargout > 1;
    for k = 1:ncols,
      if any(lindex(:,k))
        vals{k} = [find(lindex(:,k)) x(lindex(:,k),k)];
      else
        vals{k} = [];
      end
    end
    
    varargout(1) = {vals};
  end
  
