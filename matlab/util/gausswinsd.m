function window = gausswinsd(sd, Fs, ndevs)
% GAUSSWINSD Create a gaussian window with specified std. dev
%
% $Id: gausswinsd.m 2213 2009-08-03 19:38:21Z tjd $

% sd = desired standard deviation, in seconds
% Fs = sampling rate (default=1)
% ndevs = length of kernel to return, in # of stdevs (default=4)

if nargin<2
    help(mfilename)
    return
end

if nargin<3
  ndevs = 4;
end

if nargin<2
    Fs = 1;
end

npoints = round( ndevs.*sd.*Fs );
window = normpdf( linspace( -ndevs*sd, ndevs*sd, 2*npoints+1 )', 0, sd );

window = window ./ sum(window);
