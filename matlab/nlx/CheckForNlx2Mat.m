function CheckForNlx2Mat(type)
% CHECKFORNLX2MAT test that user has Neuralynx import utils installed

if nargin<1
  type = 'csc';
end

switch lower(type)
 case 'csc'
  found = exist('Nlx2MatCSC')==3;
 case 'spike'
  found = exist('Nlx2MatSpike')==3;
 case 'ev'
  found = exist('Nlx2MatEV')==3;
 case 'ts'
  found = exist('Nlx2MatTS')==3;
 case 'vt'
  found = exist('Nlx2MatVT')==3;
 otherwise
  error('Nlx2Mat filetype not recognized');
end

if ~found,  
  error(['You must install the Matlab import/export utilities from ' ...
         'Neuralynx. For PC: www.neuralynx.com/Software ; for Linux: ' ...
         'www.urut.ch/downloads/Nlx2Matv3_Aug09.tar.gz ; make sure ' ...
         'they are on your matlab path.']);
end
