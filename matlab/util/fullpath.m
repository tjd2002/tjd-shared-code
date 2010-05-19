function p = fullpath(rel_path)
%FULLPATH get full path for a given relative path
%
%  Syntax
%
%      p = fullpath( path )
%
% Does NOT search matlab's path.

if ~ischar(rel_path)
  error('Path argument must be a string');
end

if nargin<1 || isempty(rel_path)
  p = [];
  return,
end

old_path = pwd;
try
  if isdir(rel_path) % dir
    cd(rel_path);
    p = pwd;
    
  elseif ~isempty(dir(rel_path)), %file
    [pathstr, name, ext, versn] = fileparts(rel_path);
    if ~isempty(pathstr)
      cd(pathstr);
    end
    p = fullfile(pwd, [name ext versn]);
  
  else % not a dir or a file
    error(['Not a valid directory or file: ' rel_path]);
  end
  
catch
  cd(old_path);
  error(lasterr);
end

cd(old_path);
