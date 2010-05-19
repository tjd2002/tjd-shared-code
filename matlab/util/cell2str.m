function s = cell2str( c )
%CELL2STR convert cell to string representation
%
%  string=CELL2STR(cell) converts a cell array to a string
%  representation. The cell array can be recreated by eval(string). This
%  function can handle 2d cell arrys, 2d numeric and logical arrays,
%  character strings (not arrays) and structures.
%
%  See also: MAT2STR, STRUCT2STR
%

%  Copyright 2005-2008 Fabian Kloosterman

%check input arguments
if nargin<1
  help(mfilename)
  return
end

if ~iscell(c)
  error('cell2str:invalidArgument', 'Not a cell')
end

if ndims(c)>2
  error('cell2str:TwoDInput', 'Cell matrix must be 2-D.')
end

[nrows, ncols] = size( c );

%handle special cases
if isempty(c)
  if nrows==0 && ncols==0
    s = '{}';
  else
    s = ['cell(' int2str(nrows) ',' int2str(ncols) ')'];
  end
  return
end

%convert all cells
s = '{';

for k=1:nrows
  for j = 1:ncols
    if isnumeric(c{k,j}) || islogical(c{k,j})
      s = [s mat2str(c{k,j})];
    elseif iscell(c{k,j})
      s = [s cell2str(c{k,j})];
    elseif ischar(c{k,j}) && size(c{k,j},1)==1
      s = [s '''' c{k,j} '''' ];
    elseif isstruct(c{k,j})
      s = [s struct2str(c{k,j})];
    else %i.e. >1-D char arrays, function handles, objects, etc
      error('cell2str:invalidElement', ['Unable to convert ' ...
                          'object of class ' class( c{k,j} ) ' and size [' ...
                          num2str( size( c{k,j} ) ) '].'])
    end
    if j<ncols
      s = [s ','];
    end
  end
  if k<nrows
    s = [s ';'];
  end
end

s = [s '}'];
