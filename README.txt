Welcome to Tom's shared code repository!

You can either manually add the 'matlab' subdirectory and all children to 
your matlab path using the 'Path' menu, or you can add the following lines to 
your matlab startup file (e.g. ~/.matlab/startup.m). (Be sure to replace
'/path/to/thisdir' with the path to your checked-out copy of this
code.)

% Add Tom's shared Matlab code directories to Matlab's path
cd /path/to/thisdir/matlab
path(td_sharedcode_pathdef,path);


*** Directories 

matlab/
|-nlx-import
|  Code to import Neuralynx Cheetah data files into matlab structures
|
|-mwl-import
|  Code for importing Wilson Lab (MIT) data
|
|-import/
|  Code to import Wilson Lab experimental data into matlab structures.
|
|-util/
|  General-purpose utilities
|
|-cont/
|  Library of signal processing code (cdats, etc).
|
|-tdt/
|  Code to control TDT OpenEx from Matlab, and to import TDT data.
|
|-photometry/
|  Code to process recorded fiber photometry data (relies on signal 
|  processing lib in 'cont')


*** Dependencies:

 -I am currently running Matlab R2014b, and I don't try to maintain 
  backwards compatibility, so no guarantees about supporting older versions.
  (That said, code that I haven't run in a while may not be updated to work
  on latest Matlab!)

 -Matlab toolboxes:
   Filter Design Toolbox
   Image Processing Toolbox
   Signal Processing Toolbox
   Statistics Toolbox

 -mwl-import code requires Fabian Kloosterman's mwlIO library for reading and creating
  AD-style data files. (known to work with version in SVN repo as of Aug. 2009)


Any questions, or if you'd like to check in bug fixes/improvements, email
me at tjd@stanford.edu (or tjd@alum.mit.edu if you are reading this far in
the future).

Tom
Jan, 2011

