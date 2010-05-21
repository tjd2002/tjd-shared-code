Welcome to Tom's shared code repository!

You can either manually add the various subdirectories to your matlab
path using the 'Path' menu, or you can add the following lines to your
matlab startup file (e.g. ~/.matlab/startup.m). Be sure to replace
'/path/to/thisdir' with the path to your checked-out copy of this
code.

% Add Tom's shared Matlab code directories to Matlab's path

cd /path/to/thisdir/ 
path(tjd_shared_matlab_pathdef,path);

*** Directories 

matlab/
|-import/
|  Code to import experimental data into matlab structures.
|  Many dependencies on Fabian's mwlIO library
|
|-util/
|  General-purpose utilities
|
|-cont/
|  Library of signal processing code.


*** Dependencies:

 -Fabian Kloosterman's mwlIO library for reading and creating AD-style
  data files. (known to work with version in SVN repo as of Aug. 2009)

 -Matlab toolboxes. Known to work with the following versions from
  R2007a, and so far everything seems to be working with R2010a, too.
   Filter Design Toolbox >=4.1
   Image Processing Toolbox >=5.4
   Signal Processing Toolbox >=6.7
   Spline Toolbox  >=3.3.2
   Statistics Toolbox >=6.0


Any questions, or if you'd like to check in bug fixes/improvements, email
me at tjd@stanford.edu (or tjd@alum.mit.edu if you are reading this far in
the future).

Tom
May 19, 2010

