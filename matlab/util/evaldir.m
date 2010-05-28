function varargout = evaldir(epath, expr)
% EVALDIR like 'eval', but in a specified dir
%
% usage: 
% 
%   out = evaldir('/path/to/mfile.m');
%   out = evaldir('/path/to', 'mfile'); <- note no '.m'

  error(nargchk(1,2,nargin));
  
  if ~exist('expr', 'var'),
    expr = [];
  end
  
  if ~ischar(epath) && ~isempty(epath)
    error('''epath'' arg must be  string if provided');
  end

  if  ~ischar(expr) && ~isempty(expr)
    error('''expr'' arg must be  string if provided');
  end

  % Matlab likes to cache files not on the path. Clear functions fixes this 
  dbs = dbstatus;  % (keep debug state across 'clear functions')
  clear functions;
  dbstop(dbs);
  
  switch exist(epath)
   case 7 % valid directory: cd and evaluate expr
    pushd (epath);
    try
      varargout{:} = eval(expr);
    catch
      popd;
      rethrow(lasterror);
    end
    popd;

   case 2 % m-file or ordinary file on path
    [p n e] = fileparts(epath);
    if isempty(p),
      % file on matlab's path (including in current dir)
      varargout{:} = eval(epath);
    else
      % file in a specified directory
      pushd(p);
      try
        % don't include extension ('.m', etc)
        varargout{:} = eval(n);
      catch
        popd;
        error(lasterror);
      end
      popd;
    end
   
   case 0
    if isempty(epath)
      varargout{:} = eval(expr);
    else
      error(['path does not exist: ' epath]);
    end
    
   otherwise
    error('evaldir requires valid path or file for first arg');
  end